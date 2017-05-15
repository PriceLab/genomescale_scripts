#/bin/bash

# mount ephemeral drive
#sudo -s

#sudo mkfs.ext4 /dev/xvdba
#sudo mkdir -m 000 /scratch 
#echo "/dev/xvdba /scratch auto noatime 0 0" | sudo tee -a /etc/fstab
#sudo mount /scratch

# istall latest version of R
#echo "deb http://ftp.osuosl.org/pub/cran/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
echo "deb http://ftp.osuosl.org/pub/cran/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#echo "deb http://apt.postgresql.org/pub/repos/apt/ utopic-pgdg main" >> /etc/apt/sources.list
#wget -q -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -

#apt-get update && apt-get -y upgrade && apt-get -y install r-base postgresql postgresql-contrib awscli libcurl4-openssl-dev libxml2-dev libpq-dev libmariadbclient-dev && apt-get -qq install git
apt-get update && apt-get -y upgrade && apt-get -y install r-base postgresql postgresql-contrib awscli libcurl4-openssl-dev libxml2-dev libpq-dev libmariadb-client-lgpl-dev && apt-get -qq install git

Rscript setup.R

mkdir -p /scratch/data
cd /scratch/data
#aws s3 cp s3://cory-temp/gtex.primary2.rds .
#aws s3 cp s3://cory-temp/gtex.fib.RData .
#aws s3 cp s3://cory-temp/first100.RDS .
#mkdir -p /scratch/data/footprints
#cd /scratch/data/footprints
#aws s3 cp s3://cory-dbtest/footprints . --recursive

cd /scratch
mkdir -p /scratch/github
cd /scratch/github/
git clone https://github.com/PriceLab/TReNA.git
cd /scratch/github/TReNA/
git checkout 2b748798c5363c37f5a4db948c8078e1b60a34a4
# change the R version requirement
#sed -ie 's/3.4.0/3.3.0/g' /scratch/github/TReNA/DESCRIPTION
R CMD INSTALL .
aws s3 cp s3://cory-temp/coryGenomeScaleModel.R /scratch/github/TReNA/inst/utils/

mkdir -p /scratch/db
cd /scratch/db
#aws s3 cp s3://cory-dbtest/skin_hint.dump .
aws s3 cp s3://cory-dbtest/hg38.dump .
#aws s3 cp s3://cory-dbtest/fimo.dump .
aws s3 cp s3://cory-dbtest/brain_hint_16.dump .
aws s3 cp s3://cory-dbtest/brain_hint_20.dump .
aws s3 cp s3://cory-dbtest/brain_wellington_16.dump .
aws s3 cp s3://cory-dbtest/brain_wellington_20.dump .

# change the default location for postgres database storage
# stoping postgres and changing the default directory to the scratch
# stop postgres
sudo /etc/init.d/postgresql stop
mkdir -p /scratch/post
rsync -av /var/lib/postgresql/ /scratch/post/
# change default directory to /scratch/post/9.3/main/
#sed -ie 's/var\/lib\/postgresql/scratch\/post/g' /etc/postgresql/9.3/main/postgresql.conf
sed -ie 's/var\/lib\/postgresql/scratch\/post/g' /etc/postgresql/9.5/main/postgresql.conf

# increase memory size and number of connections
sed -ie 's/128MB/256MB/g' /etc/postgresql/9.5/main/postgresql.conf
sed -ie 's/max_connections\ =\ 100/max_connections\ =\ 300/g' /etc/postgresql/9.5/main/postgresql.conf
sudo /etc/init.d/postgresql start

# set up postgres
# postgres setup commands
cd /scratch/db
sudo -u postgres psql postgres << EOF
CREATE ROLE root WITH SUPERUSER CREATEDB CREATEROLE LOGIN ENCRYPTED PASSWORD 'trena';
CREATE DATABASE brain_hint_16;
CREATE DATABASE brain_hint_20;
CREATE DATABASE brain_wellington_16;
CREATE DATABASE brain_wellington_20;
CREATE DATABASE hg38;
EOF

sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=hg38 --create hg38.dump &
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=brain_hint_16 --create brain_hint_16.dump &
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=brain_hint_20 --create brain_hint_20.dump &
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=brain_wellington_16 --create brain_wellington_16.dump &
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=brain_wellington_20 --create brain_wellington_20.dump &
#sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=fimo --create fimo.dump &

wait

cd /scratch/db
sudo -u postgres psql postgres << EOF
CREATE ROLE trena WITH SUPERUSER CREATEDB CREATEROLE LOGIN ENCRYPTED PASSWORD 'trena';
grant all privileges on database brain_hint_16 to trena;
grant all privileges on database brain_hint_20 to trena;
grant all privileges on database brain_wellington_16 to trena;
grant all privileges on database brain_wellington_20 to trena;
grant all privileges on database hg38 to trena;
#grant all privileges on database fimo to trena;
EOF
