# genomescale_scripts
scripts I use for making genome-scale transcriptional regulatory networks
## This is for an i3.4xlarge instance

# mount ephemeral drive
sudo -s

sudo mkfs.ext4 /dev/xvdba
sudo mkdir -m 000 /scratch 
echo "/dev/xvdba /scratch auto noatime 0 0" | sudo tee -a /etc/fstab
sudo mount /scratch

# istall latest version of R
echo "deb http://ftp.osuosl.org/pub/cran/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
#echo "deb http://apt.postgresql.org/pub/repos/apt/ utopic-pgdg main" >> /etc/apt/sources.list
#wget -q -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -

apt-get update && apt-get -y upgrade && apt-get -y install r-base postgresql postgresql-contrib awscli libcurl4-openssl-dev libxml2-dev libpq-dev libmariadbclient-dev && apt-get -qq install git

# Set up R environment
R
source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite(c("glmnet", "GenomicRanges", "RSQLite", "lassopv", "randomForest", "flare", "vbsr", "foreach", "doParallel", "RPostgreSQL", "RMySQL", "DBI"), ask=FALSE)
q()
n

# clone TReNA
cd /scratch
mkdir -p /scratch/github
cd /scratch/github/
git clone https://github.com/PriceLab/TReNA.git
cd /scratch/github/TReNA/
R CMD INSTALL .

# specify s3 region 
# this isn't right because I have s3 privleges but the region isn't specified
aws configure
us-west-2
# get expression data
mkdir -p /scratch/data
cd /scratch/data
aws s3 cp s3://cory-temp/gtex.pc.primary.rds .

# get the sql databases
mkdir -p /scratch/db
cd /scratch/db
aws s3 cp s3://cory-dbtest/skin_hint.dump .
aws s3 cp s3://cory-dbtest/hg38.dump .

# change the default location for postgres database storage
# stoping postgres and changing the default directory to the scratch
# stop postgres
sudo /etc/init.d/postgresql stop
mkdir -p /scratch/post
rsync -av /var/lib/postgresql/ /scratch/post/
# change default directory to /scratch/post/9.3/main/
vi /etc/postgresql/9.3/main/postgresql.conf
sudo /etc/init.d/postgresql start

# set up postgres
# postgres setup commands
sudo -u postgres psql postgres
CREATE ROLE root WITH SUPERUSER CREATEDB CREATEROLE LOGIN ENCRYPTED PASSWORD 'trena';
CREATE DATABASE skin_hint;
CREATE DATABASE hg38;
\q
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=skin_hint --create skin_hint.dump &
sudo pg_restore --verbose --clean --no-acl --no-owner --dbname=hg38 --create hg38.dump &
wait
sudo -u postgres psql postgres
CREATE ROLE trena WITH SUPERUSER CREATEDB CREATEROLE LOGIN ENCRYPTED PASSWORD 'trena';
grant all privileges on database skin_hint to trena;
grant all privileges on database hg38 to trena;
\q

# set up data in R
nice R
library(TReNA)
library(doParallel)
library(RPostgreSQL)
source("/scratch/github/TReNA/inst/utils/createGenomeScaleModel.R")
gtex.primary <- readRDS("/scratch/data/gtex.primary2.rds")
gene.list <- rownames(gtex.primary)
driver <- PostgreSQL()
host <- "localhost"
dbname <- "hg38"
user <- "trena"
password <- "trena"
genome.db <- dbConnect(driver, host=host, dbname=dbname, user=user, password=password)
all.protein_coding.genes <- dbGetQuery(genome.db, "select gene_name from gtf where gene_biotype = 'protein_coding' and moleculetype = 'gene'")
pc_gene.list <- gene.list[which(gene.list %in% all.protein_coding.genes$gene_name)]
gtex.pc.primary <- gtex.primary[which(rownames(gtex.primary) %in% pc_gene.list),]









