function num_points(X, q)
    /*
    Given X / Fq returns #X(Fq) 
    //  #X(Fq) / (q+1)
    // see drews Magma ? 
    // Is input X defined over Q or Fq?
    */

    
end function;
function pwr_series_at_differentials(X, d)
    // Compute pwr series expansions of 1-forms w1..wg at degree <= d points Q1..Qn up to precision d / deg(Qi)
    // Qi: deg 1 (<= d) points in X(Fq)
    // 

end function;
function divisor_candidates(X, d)
    // generate the divisors to loop over 
    // represent divisors as pair (pt, deg)
end function;
function get_matrix(X, D)
    // use series expansions and divisor candidates to compute kernel 
    // X: curve
    // D: divisor 
end function;