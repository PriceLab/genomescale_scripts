# genomescale_scripts
scripts I use for making genome-scale transcriptional regulatory networks
## This is for an i3.4xlarge instance

# mount ephemeral drive
sudo -s

sudo mkfs.ext4 /dev/xvdba
sudo mkdir -m 000 /scratch 
echo "/dev/xvdba /scratch auto noatime 0 0" | sudo tee -a /etc/fstab
sudo mount /scratch

# enable R graphs to work through X11
# on local machine do the following
# this enables outside connections
/usr/X11R6/bin/xhost +

# when ssh into ec2 instance
ssh -XA ubuntu@amazon-instance






