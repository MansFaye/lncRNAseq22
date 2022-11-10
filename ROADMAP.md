# Roadmap
Normaly the roadmap is used to set the mission. 
I will use it also to report in more details, what I have done.

## Copy files in the projectfolder
1. get access to the cluster (modified /.ssh/config)
	> ssh idefix


2. go to the directories
Working Directories:
- /data/courses/rnaseq_course/lncRNAs/Project1/users/
- /data/courses/rnaseq_course/lncRNAs/Project2/users/

Raw Data:
- /data/courses/rnaseq/lncRNAs/fastq/

Reference files - human genome (please coordinate to download once):
- /data/courses/rnaseq/lncRNAs/Project#/references


3. Change the node (never use login node for jobs)
	> srun --cpus-per-task=1  --mem-per-cpu=500M --time=00:05:00 --pty bash


4. make my own directory for the project
	> mkdir /data/courses/rnaseq_course/lncRNAs/Project2/users/pamrein
	

5. because I realized, that everyone has the permission to change my stuff, I change them.
With Gilan I proofed, that she can't remove my folder and I can't remove her folder.
	> chmod g-w-r pamrein

HINT: Files, who are hidden are written with a dot befor the name. (example: .ssh)
To see them, you need the command ls -a


## 2022/11/06 - SSH & GIT
I set up GIT from Uni Bern and get a proper ssh connection with the ssh-agent to work.


1. Set up a new key for the ibu Cluster with the ED25519 SSH keys (book "Practical Cryptography With Go")
ssh-keygen -t ed25519 -C "pascal.amrein@students.unibe.ch"


2. copy the public key to the GIT-Lab
https://docs.gitlab.com/ee/user/ssh.html#add-an-ssh-key-to-your-gitlab-account
> cat ~/.ssh/id_ed25519.pub


3. Check if the connection works
> ssh -T git@gitlab.bioinformatics.unibe.ch

I had to give my passphrase, what means, that the ssh-agent isn't running.
I had to add the new key there. (This you have to repeat by every startup).
> eval "$(ssh-agent -s)"
> ssh-add /home/pamrein/.ssh/id_*

4. Push my files to GIT
> git add *
> git commit * -m "message"
> git push
> git pull origin master
> git push origin master 


## 08/11/2022 - Quality check of the fasta files

