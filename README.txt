
=============================================================
REST_App.py 기본 설정


MAX_RUN = 2  --> 최대 2개의 계산을 수행한다.
CPU = 55  --> Workstation에서는 큰 숫자를 하지만, 일반 서버에서는 5~10 사이의 값으로 줄인다.
TOP = 30000 --> 결과 유전자를 몇 개만 보여줄지 결정한다. 지금은 모든게 필요하게 큰 수로 잡았지만,
				서비스할때는 500이나 1000으로 제한한다.
ITERATION = 1000 -> 각 Omics / meta-data를 이용해서 p-value를 계산할 때, 몇 번의 랜덤을 이용할 것인가 결정.
					서비스할때는 100으로 줄인다.

* 주의: CPU=55로 실행하고, iteration=1000으로 하면 아래의 에러가 뜬다.
OSError: [Errno 24] Too many open files
이는 동일 파일이 너무 많이 열렸다는 뜻임.

ulimit -a를 입력하면 중간에 아래와 같은 내용이 나옴. 파일은 1024번만 열수있다는 뜻임.
open files                      (-n) 1024

이를 늘려주려면 아래의 명령을 입력. 단 세션에만 영향을 주므로 세션이 닫히면 해당 설정은 사라짐.
ulimit -n 66535

=============================================================
REST_flask_server.py

다른 건 건드릴 게 없고,
MAGIC_IP = "192.168.*.*" 만 필요시 수정하면 됨. 여기에 해당하는 IP는 최대한의 작업을 요청할 수 있음.
			그렇지 않은 경우에는 한번에 1가지 작업만 요청 가능.


pip install Werkzeug==0.16
pip install flask_restplus

python -m pip install -U py-cpuinfo

이렇게 해야 오류가 안 남.


=============================================================














Usage> 
(Multiple diseases)
python _100_NetworkExpansionModel.py -config=config_human_test.txt -cpu=10 


Output 파일이름 뒤에 _final.txt라고 붙는게 최종 결과임. 
예) TEST_208423.txt_gene.txt_final.txt


Single diseases)
python _200_NetworkExpansionModel_single_disease.py -config="./config_human_test.txt" -mod=./modifiers/human/AD.txt -cpu=10 -iteration=5 -group=2
# group: modifier를 group 수 만큼 쪼개서 계산한다.





This script runs everything!
For settings, please see config.txt


To calculate AUC, please each script with an argument 'test'. For example, to calculate GO AUC
run "python _1_GO.py test"
Output will be saved in TEMP_FOLDER (see config.txt)


[Prerequisites]
- Blast
	It can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
	I used 2.2.29.
	You should change several path informations in Blast.py

	see LOCALBLAST class
		you should change bin_path, db_path, and ofile
		and TEMP_PATH in runBLASTUsingSequence()
		
	If installed, please run Blast.py to test after chaning several settings in __main__

  
