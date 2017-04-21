#/bin/bash

R -f chunk.R 1 2170 fib.1.rds
aws s3 cp fib.1.rds s3://cory-temp/skin/

R -f chunk.R 2171 4340 fib.2.rds
aws s3 cp fib.2.rds s3://cory-temp/skin/

R -f chunk.R 4341 6510 fib.3.rds
aws s3 cp fib.3.rds s3://cory-temp/skin/

R -f chunk.R 6511 8680 fib.4.rds
aws s3 cp fib.4.rds s3://cory-temp/skin/

R -f chunk.R 8681 10850 fib.5.rds
aws s3 cp fib.5.rds s3://cory-temp/skin/

R -f chunk.R 10850 12979 fib.6.rds
aws s3 cp fib.6.rds s3://cory-temp/skin/
