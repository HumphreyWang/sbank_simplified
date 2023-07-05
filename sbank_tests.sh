export OMP_NUM_THREADS=1

# example 1 （74 templates）
python3 sbank.py --coord-frame Cartesian \
                 --x1-min 0. \
                 --x1-max 1. \
                 --distance-max 0.1 \
                 --neighborhood-size 0.2 \
                 --output-filename test_cartisian.npy \
                 --verbose \
#                 --generate-full-plots

# example 1* (22 templates)
python3 sbank.py --coord-frame Cartesian \
                 --x1-min 0. \
                 --x1-max 1. \
                 --distance-max 0.2 \
                 --neighborhood-size 0.25 \
                 --output-filename test_22.npy \
                 --verbose \

# example 2 (41 templates) (36 templates if you override the metric using [[1/4, 1/4], [1/4, 1]])
python3 sbank.py --coord-frame ScaledEuclidean \
                 --x1-min 0. \
                 --x1-max 1. \
                 --distance-max 0.1 \
                 --neighborhood-param x2 \
                 --neighborhood-size 0.2 \
                 --output-filename test_scaled.npy \
                 --verbose \
#                 --bank-seed test_22.npy \
#                 --generate-full-plots

# example 3 (54 templates)
python3 sbank.py --coord-frame Polar \
                 --x1-min 0. \
                 --x1-max 1. \
                 --x2-min 0. \
                 --x2-max 1.5707963267948966 \
                 --distance-max 0.1 \
                 --neighborhood-size 0.2 \
                 --output-filename test_polar.npy \
                 --verbose \
