function [im_obj] = make_im_obj(im, filled, iso, corners)
im_obj = struct('im', im, 'filled', filled, 'isolated', iso, 'corners', corners);
end