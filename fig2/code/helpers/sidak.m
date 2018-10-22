function alphaSID=sidak(alpha0,nc)
% Computer alpha for sidak correction, (multiple comparisons)
% alpha0= original intended alpha level
% nc = number of comparisons ("multiple comparison")
% https://en.wikipedia.org/wiki/%C5%A0id%C3%A1k_correction
% BI
alphaSID=1-((1-alpha0).^(1/nc));