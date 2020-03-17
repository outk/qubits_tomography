How to use Qubits Tomography(English)
==

I made Qubits Tomography as my graduate product, so I tell you how to use it.
src is [here](https://github.com/outk/graduate_thesis).
the theory of Qubit Tomography is written at my graduate thesis in Japanese.
So if you understand it, I recommend you read reference papers.
And I alse made Qutrit Tomography, which is three dimentional version of qubit.
If you would like to use Qudit Tomography, please use this product.
Additionally,if you can't understand these, please give it a thumbs up and comment!

As a result of this, you can get 3D graphics like this.
![mixed.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/458877/806b4508-bbcc-7f0b-3852-c8ded42b8611.png)


How to use
--

At first, copy the file of "qubit_tomo.py" in current directory.  
The structure of directory may be like this.


    .
    └project
        └qubit_tomo.py

In command prompt,

    ~\project>python qubit_tomo.py

enter "python qubit_tomo.py" or "python ./qubit_tomo.py".  

    ------------------------------------------------------------
    PLEASE ENTER NUMBER OF QUBITS
    ------------------------------------------------------------
    >>

After that, like this lines are shown, please enter a number of Qubits you want to do.  
In this time, I try 4 qubits.  
And the bases of measurement which I use is written in my graduate thesis.  
If you want to use other bases, please rewrite bases in src directly.  

     ------------------------------------------------------------
    PLEASE ENTER NUMBER OF QUBITS
    ------------------------------------------------------------
    >>
    4
    ------------------------------------------------------------
    PLEASE ENTER PATH OF EXPERIMENTAL DATA DIRECTORY

    LIKE THIS >> .\datadirectory
    ------------------------------------------------------------
    >> 

Next, please enter __directory path__ which include experimental data files ans '*.txt'.  
If there are some experimental data files, this simulater can do all of them automatically.  
But, while 3D graphics as a result is shown, simulater stop.  
So please delete graphics to continue.  
When experimental data files is like this, you can enter like this.  

    .
    └project
        ｜
        ├testdata
        ｜  ｜
        ｜  └test.txt
        ｜
        └qubit_tomo.py


    ------------------------------------------------------------
    PLEASE ENTER PATH OF EXPERIMENTAL DATA DIRECTORY

    LIKE THIS >> .\datadirectory
    ------------------------------------------------------------
    >>
    ./testdata

Next, please enter a name of result directory.  
Any name can be acceptable if your machine allow.   
If you enter nothing, the name of result directory is "default".  
I try "test".  

    ------------------------------------------------------------
    PLEASE ENTER NAME OF RESULT DIRECTORY

    THE RESULT DATA WILL SAVED AT
    '.\result\qubit\iterative(or poisson)\{ YOUR ENTED DIRECTORY NAME }\{ EXPERIMENTAL DATA FILE NAME }_result.txt'

    IF EMPTY, THE NAME OF RESULT DIRECTORY IS 'default'
    ------------------------------------------------------------
    >>
    test

After that, you need to select whether you make artificial experimental data sets or not.  
When you enter only "yes", artificial experimental data sets is made and simulated.  
These data sets is made following poisson distribution.    


    ------------------------------------------------------------
    PLEASE ENTER ANSWER WHETHER DO POISSON DISTRIBUTED SIMULATION
    IF YOU DO, PLEASE ENTER 'yes'
    IF YOU ENTER ANOTHER WORD OR EMPTY, YOUR ANSWER IS REGARED AS 'no'
    ------------------------------------------------------------
    >>
    yes
    YOUR ANSWER IS: 'yes'
    ------------------------------------------------------------
    PLEASE ENTER PATHS OF EXPERIMENTAL DATA

    IF THERE ARE MULTIPLE DATA FILE YOU WANT TO TOMOGRAPHY,
    ENTER ALL PATHS SEPARATED WITH SPACE.
    LIKE THIS >> .\datadirectory\ex1.txt .\datadirectory\ex2.txt ...
    ------------------------------------------------------------
    >>

If you enter "yes", you need to enter __a file path__ of experimental data as model.  
If you enter another word, your next question will be how many do parallel computing.  
After "yes", for example,  

    YOUR ANSWER IS: 'yes'
    ------------------------------------------------------------
    PLEASE ENTER PATHS OF EXPERIMENTAL DATA

    IF THERE ARE MULTIPLE DATA FILE YOU WANT TO TOMOGRAPHY,
    ENTER ALL PATHS SEPARATED WITH SPACE.
    LIKE THIS >> .\datadirectory\ex1.txt .\datadirectory\ex2.txt ...
    ------------------------------------------------------------
    >>
    ./testdata/test.txt

you are asked how many patterns do you make, so please enter a number you want.   

    ------------------------------------------------------------
    PLEASE ENTER ITERATION TIME OF EACH POISSON SIMULATION
    ------------------------------------------------------------
    >>
    5

Finally, you need to enter how many parallel computing you want.  
And the maximum number of parallel that your computer can do is shown, so please refer it.  
(__※Attention:__　If you choose too many number of parallel, CPU can't work well. So I recommend about half of that number. And you must check your CPU working rate is proper.)  

    ------------------------------------------------------------
    HOW MANY TIMES DO YOU WANT TO PARALLELIZE?
    IF THE NUMBER IS TOO LARGE, THE PARFORMANCE OF SIMULATION BECOME LOWER.
    THE NUMBER OF LOGICAL PROCESSOR OF YOUR COMPUTER IS >>
    6
    RECOMENDED NUMBER IS LESS THAN THE ABOVE NUMBER.
    ------------------------------------------------------------
    >>
    2

All that is left is waiting for finishing simulating.  
As mentioned above, whem 3D graphics is shown, simulater is stop.  
So if you do not need to show this picture, please comment out of line 310 in src.  


    plotResult(numberOfQubits, estimatedDensityMatrix, baseNames)

Like this.

    #plotResult(numberOfQubits, estimatedDensityMatrix, baseNames)


Result
--

The result is saved in result directory of which you enter the name.  

    .
    └project
        ｜
        ├result
        ｜  ｜
        ｜  └qubit
        ｜      ｜
        ｜      ├iterative
        ｜      ｜  ｜
        ｜      ｜  └test
        ｜      ｜      ｜
        ｜      ｜      └result.txt
        ｜      ｜
        ｜      └poisson
        ｜          ｜
        ｜          └test
        ｜              ｜
        ｜              └result.txt
        ｜
        ├testdata
        ｜  ｜
        ｜  └test.txt
        ｜
        └qubit_tomo.py

In `result.txt`, the fidelity between result density matrix and ideal matrix which you prepared.  
Ideal matrix is prepared in src directly.  
So if you need, please fix your own.  



Qubitsトモグラフィーの使い方(日本語)
==

まず、qubit_tomo.pyを作業Directory(project)下にコピーします。


    .
    └project
        └qubit_tomo.py

そこで

    ~\project>python qubit_tomo.py

と入力します。  

    ------------------------------------------------------------
    PLEASE ENTER NUMBER OF QUBITS
    ------------------------------------------------------------
    >>

と表示されるのでシミュレートしたいQubitsの数を入力してください。  
今回は4つで試します。  
また、測定基底は卒業論文に記載しているとおりです。  
他の基底で測定した場合はcodeを書き換えてください。  

     ------------------------------------------------------------
    PLEASE ENTER NUMBER OF QUBITS
    ------------------------------------------------------------
    >>
    4
    ------------------------------------------------------------
    PLEASE ENTER PATH OF EXPERIMENTAL DATA DIRECTORY

    LIKE THIS >> .\datadirectory
    ------------------------------------------------------------
    >> 

次にtxtファイルで保存した実験データがあるディレクトリパスを入力してください。  
実験データが複数あってもすべてについてトモグラフィーしてくれます。  
ただし、3D描画画面が出たままではシミュレートは進みません。  
注意してください。
実験データを次のように保存している場合は以下のように指定できます。

    .
    └project
        ｜
        ├testdata
        ｜  ｜
        ｜  └test.txt
        ｜
        └qubit_tomo.py


    ------------------------------------------------------------
    PLEASE ENTER PATH OF EXPERIMENTAL DATA DIRECTORY

    LIKE THIS >> .\datadirectory
    ------------------------------------------------------------
    >>
    ./testdata

次に、計算結果の出力先ディレクトリ名を入力してください。  
ここはなんでもいいです。  
何も入力しなければ”default”になります。  
例えば、

    ------------------------------------------------------------
    PLEASE ENTER NAME OF RESULT DIRECTORY

    THE RESULT DATA WILL SAVED AT
    '.\result\qubit\iterative(or poisson)\{ YOUR ENTED DIRECTORY NAME }\{ EXPERIMENTAL DATA FILE NAME }_result.txt'

    IF EMPTY, THE NAME OF RESULT DIRECTORY IS 'default'
    ------------------------------------------------------------
    >>
    test

ここで、疑似実験データを作成するか聞かれます。  
"yes"以外の入力はすべて"no"と判断されます。  
疑似実験データは各実験データの回数をそれぞれ期待値としたポアソン分布に沿うようにランダムに生成されます。  


    ------------------------------------------------------------
    PLEASE ENTER ANSWER WHETHER DO POISSON DISTRIBUTED SIMULATION
    IF YOU DO, PLEASE ENTER 'yes'
    IF YOU ENTER ANOTHER WORD OR EMPTY, YOUR ANSWER IS REGARED AS 'no'
    ------------------------------------------------------------
    >>
    yes
    YOUR ANSWER IS: 'yes'
    ------------------------------------------------------------
    PLEASE ENTER PATHS OF EXPERIMENTAL DATA

    IF THERE ARE MULTIPLE DATA FILE YOU WANT TO TOMOGRAPHY,
    ENTER ALL PATHS SEPARATED WITH SPACE.
    LIKE THIS >> .\datadirectory\ex1.txt .\datadirectory\ex2.txt ...
    ------------------------------------------------------------
    >>

"yes"と入力すると、このように表示されるので、疑似実験データのもととなる実験データ**ファイルパス**を入力してください。  
"no"の場合は最後のの並列化数入力まで飛ばされます。    
"yes"の場合、例えば、

    YOUR ANSWER IS: 'yes'
    ------------------------------------------------------------
    PLEASE ENTER PATHS OF EXPERIMENTAL DATA

    IF THERE ARE MULTIPLE DATA FILE YOU WANT TO TOMOGRAPHY,
    ENTER ALL PATHS SEPARATED WITH SPACE.
    LIKE THIS >> .\datadirectory\ex1.txt .\datadirectory\ex2.txt ...
    ------------------------------------------------------------
    >>
    ./testdata/test.txt

すると、疑似実験データを何パターン生成するか聞かれるので生成したい数を入力してください。  

    ------------------------------------------------------------
    PLEASE ENTER ITERATION TIME OF EACH POISSON SIMULATION
    ------------------------------------------------------------
    >>
    5

最後に、いくつ並列化させて計算するか聞かれます。  
また、現在使用しているパソコンの並列化可能な最大の数が表示されるので参考にしてください。  
（ **※注意：** 並列化可能な最大の数の半分程度でないとCPU使用率が跳ね上がってシミュレートが進まなくなります。シミュレート実行後に、必ずリソースモニターなどでCPU使用率を確認してください。）

    ------------------------------------------------------------
    HOW MANY TIMES DO YOU WANT TO PARALLELIZE?
    IF THE NUMBER IS TOO LARGE, THE PARFORMANCE OF SIMULATION BECOME LOWER.
    THE NUMBER OF LOGICAL PROCESSOR OF YOUR COMPUTER IS >>
    6
    RECOMENDED NUMBER IS LESS THAN THE ABOVE NUMBER.
    ------------------------------------------------------------
    >>
    2

あとは、シミュレートが終わるまで待つだけです。  
上でも述べましたが、3D描画が表示されている間は他のシミュレートは進まないので、必要がない場合はsrcのl.310の

    plotResult(numberOfQubits, estimatedDensityMatrix, baseNames)

を次のようにコメントアウトしてください。

    #plotResult(numberOfQubits, estimatedDensityMatrix, baseNames)


出力結果
--

出力結果は次のディレクトリに保存されます。

    .
    └project
        ｜
        ├result
        ｜  ｜
        ｜  └qubit
        ｜      ｜
        ｜      ├iterative
        ｜      ｜  ｜
        ｜      ｜  └test
        ｜      ｜      ｜
        ｜      ｜      └result.txt
        ｜      ｜
        ｜      └poisson
        ｜          ｜
        ｜          └test
        ｜              ｜
        ｜              └result.txt
        ｜
        ├testdata
        ｜  ｜
        ｜  └test.txt
        ｜
        └qubit_tomo.py

`result.txt`には理想状態とのfidelityが保存されています。  
理想状態は`qubit_tomo.py`に直接用意しています。  
必要に応じて書き換えてください。