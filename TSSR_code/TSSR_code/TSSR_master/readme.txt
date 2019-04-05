--method 			set LDA (lncRNA-disease association)prediction method
--dataset: 			choose the benchmark dataset, i.e., lncRNADisease,mndr,lnc2cancer
--folder:			set the the folder that contains the datasets (default "datasets
--specify-arg:		 1 for using default/specified arguments (default 1)
--method-opt:		set arguments for each method (method ARGUMENTS have the form name=value)
--predict-num:		a positive integer for predicting top-N novel LDAs (default 100)

some examples:

(1) run a method with default arguments
	python PyLDA.py --method="tssr" --dataset="mndr"
	python PyLDA.py --method="tssr" --dataset="mndr" --predict-num=100
	python PyLDA.py --method="tssr" --dataset="mndr"  --specify-arg=1
	python PyLDA.py --method="tssr" --dataset="mndr"  --specify-arg=1 --predict-num=100

(2) run a method with specified arguments

	python PyLDA.py --method="tssr" --dataset="mndr" --specify-arg=1 --method-opt="beta=4"

	python PyLDA.py --method="tssr" --dataset="mndr" --specify-arg=1 --method-opt=" lambda_d=1  lambda_t=1 beta=256"


(3) predict the top-N novel LDAs

	python PyLDA.py --method="tssr" --dataset="mndr" --predict-num=N --method-opt="lambda_d=1 lambda_t=1  beta=256"