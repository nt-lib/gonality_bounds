import functools
import time
from lmf import db

# Computes the inequality in the Castelnuovo-Severi theorem. If the gonality is below this bound, then the gonal map
# and the map X_K -> X_H both factor through some X_L
def castelnuovo_severi_inequality(K_genus, H_genus, degree):
    return (K_genus - degree*H_genus + degree - 2)//(degree - 1) + 1

# Computes the Castelnuovo-Severi bound for the gonality of X_K given a modular inclusion map X_K -> X_H
# Requires that K is a maximal subgroup of H
def castelnuovo_severi(K_genus, H_genus, degree, H_gonality_lower_bound):
    cs_bound = castelnuovo_severi_inequality(K_genus, H_genus, degree)
    gonality_bound = degree*H_gonality_lower_bound
    return min(cs_bound, gonality_bound)

level_bound = 40
print("Loading modular curves of level <=", level_bound, "from the LMFDB")
curves = {x["label"]: x for x in db.gps_gl2zhat_fine.search({"contains_negative_one": True, "level": {"$lte": level_bound}}, ["label", "index", "genus", "parents", "q_gonality_bounds"])}
print("Number of curves:", len(curves))

sorted_labels = sorted(curves.keys(), key=lambda x: curves[x]["index"])
for x in curves.values():
    x["gonality_lower_bound"] = x["q_gonality_bounds"][0]

while True:
    print("\n===== Applying C-S to minimal covers =====")
    start_time = time.time()
    unknown_gonality_count = 0
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in sorted_labels:
        K = curves[K_label]
        if K["gonality_lower_bound"] != K["q_gonality_bounds"][1]:
            unknown_gonality_count += 1
            max_cs_bound = 0
            for H_label in K["parents"]:
                H = curves[H_label]
                d = K["index"] // H["index"]
                cs_bound = castelnuovo_severi(K["genus"], H["genus"], d, H["gonality_lower_bound"])
                max_cs_bound = max(max_cs_bound, cs_bound)
            if max_cs_bound > K["gonality_lower_bound"]:
                improved_gonality_count += 1
                K["gonality_lower_bound"] = max_cs_bound
            if max_cs_bound >= K["q_gonality_bounds"][1]:
                exact_gonality_count += 1

    print("Inexact gonality:", unknown_gonality_count)
    print("C-S improves bound:", improved_gonality_count)
    print("C-S obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break

    print("\n===== Propagating downwards =====")
    start_time = time.time()
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in reversed(sorted_labels):
        K = curves[K_label]
        for H_label in K["parents"]:
            H = curves[H_label]
            d = K["index"] // H["index"]
            bound = (K["gonality_lower_bound"] + d - 1)//d
            if bound > H["gonality_lower_bound"]:
                improved_gonality_count += 1
                H["gonality_lower_bound"] = bound
                if bound >= K["q_gonality_bounds"][1]:
                    exact_gonality_count += 1
    print("Downwards propagation improves existing bound:", improved_gonality_count)
    print("Downwards propagation obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break

while True:
    print("\n===== Applying C-S to all covers =====")
    start_time = time.time()
    unknown_gonality_count = 0
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in sorted_labels:
        K = curves[K_label]
        if K["gonality_lower_bound"] != K["q_gonality_bounds"][1]:
            unknown_gonality_count += 1
            if len(K["parents"]) == 0:
                # No curves to use Castelnuovo-Severi
                continue

            max_cs_bound = 0
            queue = {}
            factoring_bounds = {}
            for L_label in K["parents"]:
                L = curves[L_label]
                if L["index"] not in queue:
                    queue[L["index"]] = {L_label}
                else:
                    queue[L["index"]].add(L_label)
                factoring_bounds[L_label] = K["index"]//L["index"]*L["gonality_lower_bound"]
            current_index = max(queue.keys())
            while True:
                # As Abbey pointed out, can avoid exploring nodes where the gonality bound is already known to be less
                # than the current max_cs_bound. How does one implement this pruning in a non-tree?
                if len(queue[current_index]) == 0:
                    queue.pop(current_index)
                    if len(queue) == 0:
                        break
                    current_index = max(queue.keys())
                H_label = queue[current_index].pop()
                H = curves[H_label]
                d = K["index"] // H["index"]
                cs_inequality = castelnuovo_severi_inequality(K["genus"], H["genus"], d)
                cs_bound = min(cs_inequality, factoring_bounds[H_label])
                max_cs_bound = max(max_cs_bound, cs_bound)
                
                for parent in H["parents"]:
                    index = curves[parent]["index"]
                    if index not in queue:
                        queue[index] = {parent}
                        factoring_bounds[parent] = factoring_bounds[H_label]
                    elif parent not in queue[index]:
                        queue[index].add(parent)
                        factoring_bounds[parent] = factoring_bounds[H_label]
                    elif factoring_bounds[parent] > factoring_bounds[H_label]:
                        factoring_bounds[parent] = factoring_bounds[H_label]
            if max_cs_bound > K["gonality_lower_bound"]:
                improved_gonality_count += 1
                K["gonality_lower_bound"] = max_cs_bound
            if max_cs_bound >= K["q_gonality_bounds"][1]:
                exact_gonality_count += 1

    print("Inexact gonality:", unknown_gonality_count)
    print("C-S improves bound:", improved_gonality_count)
    print("C-S obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break
        
    print("\n===== Propagating downwards =====")
    start_time = time.time()
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in reversed(sorted_labels):
        K = curves[K_label]
        for H_label in K["parents"]:
            H = curves[H_label]
            d = K["index"] // H["index"]
            bound = (K["gonality_lower_bound"] + d - 1)//d
            if bound > H["gonality_lower_bound"]:
                improved_gonality_count += 1
                H["gonality_lower_bound"] = bound
                if bound >= K["q_gonality_bounds"][1]:
                    exact_gonality_count += 1
    print("Downwards propagation improves existing bound:", improved_gonality_count)
    print("Downwards propagation obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break

while True:
    print("\n===== Applying gcd C-S to all covers =====")
    start_time = time.time()
    unknown_gonality_count = 0
    can_improve_count = 0
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in sorted_labels:
        K = curves[K_label]
        K["unconstrained_gonality_bound"] = K["gonality_lower_bound"]
        K["low_possible_gonalities"] = set()
        if K["gonality_lower_bound"] != K["q_gonality_bounds"][1]:
            if len(K["parents"]) == 0:
                # No curves to use Castelnuovo-Severi
                continue
            unknown_gonality_count += 1

            queue = {}
            for L_label in K["parents"]:
                L = curves[L_label]
                if L["index"] not in queue:
                    queue[L["index"]] = {L_label}
                else:
                    queue[L["index"]].add(L_label)
            current_index = max(queue.keys())
            while True:
                if len(queue[current_index]) == 0:
                    queue.pop(current_index)
                    if len(queue) == 0:
                        break
                    current_index = max(queue.keys())
                H_label = queue[current_index].pop()
                H = curves[H_label]
                d = K["index"] // H["index"]
                cs_inequality = castelnuovo_severi_inequality(K["genus"], H["genus"], d)
                if K["unconstrained_gonality_bound"] < cs_inequality:
                    K["unconstrained_gonality_bound"] = cs_inequality
                
                for parent in H["parents"]:
                    index = curves[parent]["index"]
                    if index not in queue:
                        queue[index] = {parent}
                    elif parent not in queue[index]:
                        queue[index].add(parent)

            if K["unconstrained_gonality_bound"] <= K["gonality_lower_bound"]:
                continue

            can_improve_count += 1
            K["low_possible_gonalities"] = set(range(K["gonality_lower_bound"], K["unconstrained_gonality_bound"]))

            queue = {}
            factoring_bounds = {}
            for L_label in K["parents"]:
                L = curves[L_label]
                if L["index"] not in queue:
                    queue[L["index"]] = {L_label}
                else:
                    queue[L["index"]].add(L_label)
                factoring_bounds[L_label] = set()
                degree = K["index"]//L["index"]
                for s in L["low_possible_gonalities"]:
                    if s*degree < K["unconstrained_gonality_bound"]:
                        factoring_bounds[L_label].add(s*degree)
                factoring_bounds[L_label].update([degree*i for i in range(L["unconstrained_gonality_bound"], (K["unconstrained_gonality_bound"] - 1)//degree + 1)])
            current_index = max(queue.keys())
            while True:
                if len(queue[current_index]) == 0:
                    queue.pop(current_index)
                    if len(queue) == 0:
                        break
                    current_index = max(queue.keys())
                H_label = queue[current_index].pop()
                H = curves[H_label]
                d = K["index"] // H["index"]
                cs_inequality = castelnuovo_severi_inequality(K["genus"], H["genus"], d)
                K["low_possible_gonalities"].intersection_update(factoring_bounds[H_label].union(range(cs_inequality, K["unconstrained_gonality_bound"])))
                
                for parent in H["parents"]:
                    index = curves[parent]["index"]
                    if index not in queue:
                        queue[index] = {parent}
                        factoring_bounds[parent] = factoring_bounds[H_label].copy()
                    elif parent not in queue[index]:
                        queue[index].add(parent)
                        factoring_bounds[parent] = factoring_bounds[H_label].copy()
                    elif factoring_bounds[parent] < factoring_bounds[H_label]:
                        factoring_bounds[parent].update(factoring_bounds[H_label])
            new_lower_bound = min(K["low_possible_gonalities"], default=K["unconstrained_gonality_bound"])
            if new_lower_bound > K["gonality_lower_bound"]:
                print(K["label"], K["low_possible_gonalities"], K["unconstrained_gonality_bound"], K["gonality_lower_bound"])
                improved_gonality_count += 1
                K["gonality_lower_bound"] = new_lower_bound
            if new_lower_bound >= K["q_gonality_bounds"][1]:
                exact_gonality_count += 1

    print("Inexact gonality:", unknown_gonality_count)
    print("C-S can obtain improvement:", can_improve_count)
    print("C-S improves bound:", improved_gonality_count)
    print("C-S obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break
        
    print("\n===== Propagating downwards =====")
    start_time = time.time()
    improved_gonality_count = 0
    exact_gonality_count = 0
    for K_label in reversed(sorted_labels):
        K = curves[K_label]
        for H_label in K["parents"]:
            H = curves[H_label]
            d = K["index"] // H["index"]
            bound = (K["gonality_lower_bound"] + d - 1)//d
            if bound > H["gonality_lower_bound"]:
                improved_gonality_count += 1
                H["gonality_lower_bound"] = bound
                if bound >= K["q_gonality_bounds"][1]:
                    exact_gonality_count += 1
    print("Downwards propagation improves existing bound:", improved_gonality_count)
    print("Downwards propagation obtains exact bound:", exact_gonality_count)
    print("Time elapsed:", time.time() - start_time)

    if improved_gonality_count == 0:
        break
