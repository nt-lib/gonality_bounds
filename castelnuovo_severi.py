from lmf import db

# Computes the Castelnuovo-Severi bound for the gonality of X_K given a modular inclusion map X_K -> X_H
# Requires that K is a maximal subgroup of H
def castelnuovo_severi(K_genus, H_genus, degree, H_gonality_lower_bound):
    cs_bound = (K_genus - degree*H_genus + degree - 2)//(degree - 1) + 1
    gonality_bound = degree*H_gonality_lower_bound
    return min(cs_bound, gonality_bound)

curves = {x["label"]: x for x in db.gps_gl2zhat_fine.search({"contains_negative_one": True, "level": {"$lte": 40}}, ["label", "index", "genus", "parents", "q_gonality_bounds"])}
print("Number of curves:", len(curves))

unknown_gonality_count = 0
improved_gonality_count = 0
exact_gonality_count = 0
for K in curves.values():
    if K["q_gonality_bounds"][0] != K["q_gonality_bounds"][1]:
        unknown_gonality_count += 1
        max_cs_bound = 0
        for H_label in K["parents"]:
            H = curves[H_label]
            d = K["index"] // H["index"]
            cs_bound = castelnuovo_severi(K["genus"], H["genus"], d, H["q_gonality_bounds"][0])
            max_cs_bound = max(max_cs_bound, cs_bound)
        if max_cs_bound > K["q_gonality_bounds"][0]:
            improved_gonality_count += 1
        if max_cs_bound >= K["q_gonality_bounds"][1]:
            exact_gonality_count += 1

print("Inexact gonality in LMFDB:", unknown_gonality_count)
print("C-S improves existing bound:", improved_gonality_count)
print("C-S obtains exact bound:", exact_gonality_count)
