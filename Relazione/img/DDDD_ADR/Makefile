file_eps := $(wildcard *.eps)
file_png := $(wildcard *.png)
file := $(subst .png,,$(file_png))
numerofile := Numbers of files png

clean:
	rm -f $(file_eps)
files:
	echo "$(numerofile) $(words $(file_png))"	
eps:
	convert $(word 1,$(file_png)) $(word 1,$(subst .png,.eps,$(file_png)))
	convert $(word 2,$(file_png)) $(word 2,$(subst .png,.eps,$(file_png)))
	convert $(word 3,$(file_png)) $(word 3,$(subst .png,.eps,$(file_png)))
	convert $(word 4,$(file_png)) $(word 4,$(subst .png,.eps,$(file_png)))
	convert $(word 5,$(file_png)) $(word 5,$(subst .png,.eps,$(file_png)))
	convert $(word 6,$(file_png)) $(word 6,$(subst .png,.eps,$(file_png)))
