Be sure that the following are installed and setup properly
	mountainlab
	nodejs
	If needed do this:
		> sudo ln -s /usr/bin/nodejs /usr/bin/node

For the tabulation step you need the console.table module. Run:
	> npm install console.table

000 - Run the processing daemon in a separate terminal
	Open a separate terminal and run:
	> mountainprocess daemon-start
	Keep that terminal open as it will handle the queing/running of scripts and processes

001 - Generate a collection of synthetic example datasets
	The raw data go into BIGFILES/synth_examples directory
	The .prv links to these go into the datasets/synth_examples directory

002 - 
