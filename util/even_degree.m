function degree = even_degree(p,vars)
  % deg = even_degree(p,vars)
  % calculate the degree of polynomial p, in terms of the variables vars
  % Rounds up to the nearest even number
  if isnumeric(p)
    degree = 0;
  else
    degree = full(2*ceil(deg(p,vars)/2));
  end
end