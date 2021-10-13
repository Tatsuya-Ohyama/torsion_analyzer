# torsion_analyzer.py

## 概要
トラジェクトリや PDB ファイルの構造内の特定構造の二面角を測定するプログラム


## 使用方法
```sh
$ torsion_analyzer.py [-h] [-p TOPOLOGY_FILE] -x COORDINATE_FILE -ma AMBERMASK_AXIS [AMBERMASK_AXIS ...] -mb AMBERMASK_TORSION [AMBERMASK_TORSION ...] -o OUTPUT_FILE [-b START_TIME] [-e END_TIME] [--offset OFFSET] [-O]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-p TOPOLOGY_FILE`
	: トポロジーファイル (.pdb for .pdb, .prmtop for .nc, and .gro for .xtc)
* `-x COORDINATE_FILE`
	: 座標 (トラジェクトリ) ファイル (.pdb, .nc and .xtc)
* `-ma AMBERMASK_AXIS [AMBERMASK_AXIS ...]`
	: 中心軸を算出するための原子の Ambermask リスト (1 原子マスク: 該当原子を中心軸候補とする / 2 原子マスク: 中心軸算出のための中心点を算出用の 2 原子 (例: `:1@CA,:10@CA`))
* `-mb AMBERMASK_TORSION [AMBERMASK_TORSION ...]`
	: 中心軸から二面角を算出するための Ambermask リスト
* `-o OUTPUT_FILE`
	: 出力ファイル (.csv or .txt (tsv))
* `-b START_TIME`
	: トラジェクトリファイルの開始フレーム (start from 0)
* `-e END_TIME`
	: トラジェクトリの終了フレーム (start from 0)
* `--offset OFFSET`
	: トラジェクトリの読み込み間隔 (Default: 1)
* `-O`
	: 上書きプロンプトを出さずに出力する。


## 動作要件
* Python3
	* numpy
	* parmed
	* MDAnalysis
	* netCDF4
	* tqdm


## License
The MIT License (MIT)

Copyright (c) 2021 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 1.0 (2021-10-13)
* 公開した。
