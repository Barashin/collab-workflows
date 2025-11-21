# ワークフロー各ノードの説明

このディレクトリには、ワークフロー画像の各ノードに対応するPythonスクリプトが含まれています。

## ノード一覧

### Node 1: PDBファイルダウンロード
- **ファイル**: `protein/node_01_download_pdb.py`
- **入力**: `PARAM_PDB_ID` (環境変数)
- **出力**: `protein/output/{pdb_id}.pdb`
- **説明**: RCSB PDBからタンパク質構造ファイルをダウンロード

### Node 2: チェーン抽出
- **ファイル**: `preparation/node_02_extract_chain.py`
- **入力**: `protein/output/{pdb_id}.pdb`
- **出力**: `preparation/output/{pdb_id}_A_NAD.pdb`
- **説明**: PDBファイルからChain AとNAD補因子を抽出

### Node 3: リガンド中心計算とconfig生成
- **ファイル**: `preparation/node_03_ligand_center_config.py`
- **入力**: `protein/output/{pdb_id}.pdb`, `PARAM_LIGAND_NAME` (環境変数)
- **出力**: `preparation/output/config.txt`
- **説明**: リガンドの中心座標を計算し、ドッキング設定ファイルを生成

### Node 4: タンパク質構造修正
- **ファイル**: `preparation/node_04_fix_structure.py`
- **入力**: `preparation/output/{pdb_id}_A_NAD.pdb`
- **出力**: `preparation/output/{pdb_id}_A_NAD_fixed_with_NAD.pdb`
- **説明**: PDBFixerを使用してタンパク質構造をクリーンアップ（欠損残基の追加、水素の追加など）し、NADを再付与

### Node 5: リガンド抽出
- **ファイル**: `ligand/node_05_extract_ligand.py`
- **入力**: `preparation/output/{pdb_id}_A_NAD.pdb`, `PARAM_LIGAND_NAME` (環境変数)
- **出力**: `ligand/output/{ligand_name}.sdf`
- **説明**: PDBファイルからリガンドを抽出してSDF形式で保存（元の座標を保持）

### Node 6: リガンドバリアント生成
- **ファイル**: `ligand/node_06_generate_variants.py`
- **入力**: `ligand/output/{ligand_name}.sdf`
- **出力**: 
  - `ligand/output/ligand_library/` ディレクトリ（複数のSDFファイル）
  - `ligand/output/variants.svg`
- **説明**: リガンドのSMILESから機能基バリアントを生成し、3D構造を生成

### Node 7: レセプター準備
- **ファイル**: `docking/node_07_prepare_receptor.py`
- **入力**: `preparation/output/{pdb_id}_A_NAD_fixed_with_NAD.pdb`
- **出力**: `docking/output/{pdb_id}_A_NAD_fixed_with_NAD.pdbqt`
- **説明**: タンパク質構造をPDBQT形式に変換（AutoDockTools使用）

### Node 8: ドッキングスクリーニング
- **ファイル**: `docking/node_08_docking_screening.py`
- **入力**: 
  - `docking/output/{pdb_id}_A_NAD_fixed_with_NAD.pdbqt`
  - `preparation/output/config.txt`
  - `ligand/output/ligand_library/` (個別のSDFファイル)
- **出力**: `docking/output/` ディレクトリ（ドッキング結果ファイル群）
- **説明**: SMINAを使用して分子ドッキングを実行

### Node 9: レポート生成
- **ファイル**: `report/node_09_reporting.py`
- **入力**: `docking/output/` ディレクトリ内のドッキング結果
- **出力**: 
  - `report/output/docking_ranking.txt`
  - `report/output/*_docked.sdf` (トップ化合物)
- **説明**: ドッキング結果を解析してランキングを生成

## 実行順序

ワークフローは以下の順序で実行する必要があります：

```
1. Node 1 (PDBダウンロード)
   ↓
2. Node 2 (チェーン抽出) ← Node 1
   ↓
3. Node 3 (リガンド中心計算とconfig生成) ← Node 1
4. Node 4 (タンパク質構造修正) ← Node 2
5. Node 5 (リガンド抽出) ← Node 2
   ↓
6. Node 6 (リガンドバリアント生成) ← Node 5
   ↓
7. Node 7 (レセプター準備) ← Node 4
   ↓
8. Node 8 (ドッキングスクリーニング) ← Node 3, 6, 7
   ↓
9. Node 9 (レポート生成) ← Node 8
```

## 実行方法

### 一括実行（推奨）

ワークフロー全体を3つのステージに分けて実行：

```bash
# ステージ1: データ準備（タンパク質のダウンロードと準備、リガンドの抽出とバリアント生成）
bash pre_run.sh

# ステージ2: ドッキングスクリーニング
bash run.sh

# ステージ3: レポート生成
bash post_run.sh
```

### 個別実行

各ノードを個別に実行する場合：

```bash
# タンパク質関連
python3 protein/node_01_download_pdb.py

# 準備関連
python3 preparation/node_02_extract_chain.py
python3 preparation/node_03_ligand_center_config.py
python3 preparation/node_04_fix_structure.py

# リガンド関連
python3 ligand/node_05_extract_ligand.py
python3 ligand/node_06_generate_variants.py

# ドッキング
python3 docking/node_07_prepare_receptor.py
python3 docking/node_08_docking_screening.py

# レポート
python3 report/node_09_reporting.py
```

## 環境変数

- `PARAM_PDB_ID`: PDB IDを指定（例: "4OHU"）
- `PARAM_LIGAND_NAME`: リガンドの3文字コードを指定（例: "2TK"）
- `PARAM_EXHAUSTIVENESS`: ドッキングの網羅性（デフォルト: "8"）
- `PARAM_NUM_MODES`: ドッキングモード数（デフォルト: "1"）
- `PARAM_ENERGY_RANGE`: エネルギー範囲（デフォルト: "3"）

```bash
export PARAM_PDB_ID="4OHU"
export PARAM_LIGAND_NAME="2TK"
export PARAM_EXHAUSTIVENESS="8"
export PARAM_NUM_MODES="1"
export PARAM_ENERGY_RANGE="3"
```

## ファイル構造

```
update_job1/
├── pre_run.sh              # データ準備ステージ
├── run.sh                  # ドッキングスクリーニングステージ
├── post_run.sh             # レポート生成ステージ
├── @job.toml               # ジョブ設定ファイル
├── protein/                # タンパク質処理ノード
│   ├── node_01_download_pdb.py
│   └── output/             # タンパク質処理の出力
├── preparation/            # 準備処理ノード
│   ├── node_02_extract_chain.py
│   ├── node_03_ligand_center_config.py
│   ├── node_04_fix_structure.py
│   └── output/             # 準備処理の出力
├── ligand/                 # リガンド処理ノード
│   ├── node_05_extract_ligand.py
│   ├── node_06_generate_variants.py
│   └── output/             # リガンド処理の出力
│       └── ligand_library/ # 生成されたリガンドバリアント
├── docking/                # ドッキングノード
│   ├── node_07_prepare_receptor.py
│   ├── node_08_docking_screening.py
│   └── output/             # ドッキング結果
└── report/                 # レポートノード
    ├── node_09_reporting.py
    └── output/             # レポート出力
```

## 注意事項

1. **Node 3 (リガンド中心計算)**: 元のPDBファイル（Node 1の出力）からリガンド座標を取得します。これは、チェーン抽出後のファイルにはリガンドが含まれていないためです。

2. **Node 4 (タンパク質構造修正)**: PDBFixerを使用して構造を修正しますが、NAD補因子が削除される可能性があるため、自動的に再付与されます。

3. **Node 6 (リガンドバリアント生成)**: SMILESから機能基バリアントを生成し、3D構造を自動生成します。元のリガンドも`ligand_library/original.sdf`として保存されます。

4. **Node 7 (レセプター準備)**: AutoDockToolsの`prepare_receptor4.py`が必要です。MGLTools/AutoDockToolsがインストールされている必要があります。

5. **Node 8 (ドッキングスクリーニング)**: SMINAがインストールされている必要があります。`ligand_library/`ディレクトリ内のすべてのSDFファイルに対してドッキングを実行します。

6. **依存関係**: 各ノードは前のノードの出力を入力として使用するため、順序を守って実行してください。

7. **出力ディレクトリ**: 
   - `protein/output/`: タンパク質処理の出力
   - `preparation/output/`: タンパク質準備の出力
   - `ligand/output/`: リガンド処理の出力
   - `docking/output/`: ドッキング結果
   - `report/output/`: 最終結果（ランキングと最良のリガンド）

