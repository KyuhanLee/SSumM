all: compile demo
compile:
	-chmod u+x ./*.sh
	./compile.sh
demo:
	-chmod u+x ./*.sh
	rm -rf output
	mkdir output
	@echo [DEMO] running SSumM
	./run.sh ego_facebook.txt 0.20 1
	@echo [DEMO] summary graph is saved in the output directory
