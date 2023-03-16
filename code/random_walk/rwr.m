function W = rwr(P, restart_prob)

  n = size(P, 1); 
  W = (eye(n) - (1 - restart_prob) * P) \ (restart_prob * eye(n));
  
end