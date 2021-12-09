function [ntest, test_filt] = cv_partition(anno, test_frac)
    ng = size(anno, 2);
    ntest = floor(ng * test_frac);
    test_ind = randperm(ng, ntest);
    test_filt = false(ng, 1);
    test_filt(test_ind) = true;
    ntest = length(test_ind);
end