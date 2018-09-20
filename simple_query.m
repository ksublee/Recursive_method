function q_index = simple_query(a, b)
    % b should be monotone increasing with equal distance
    % if query points are out of bound, take maximum or minimum
    a(a > max(b)) = max(b);
    a(a < min(b)) = min(b);
    dx = b(2) - b(1);
    q_index = round((a - min(b)) / dx) + 1;  % index
end
