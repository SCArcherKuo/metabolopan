# 使用手冊（User Manual）

本手冊記錄本軟體在數值上的運作方式——演算法、預設門檻值、以及與常見替代做法的差異——讓你能在論文或報告中為它產生的每一個數字辯護。在發表任何依賴本軟體的結果之前，請先完整閱讀一次。

> 本文件為 [`USER_MANUAL.md`](USER_MANUAL.md) 的台灣繁體中文翻譯版。若中英文版本有出入，請以英文版為準。

本軟體的分析流程分為三個階段。預設門檻值以方括號 `[...]` 標示。

- [第一階段 — 輸入解析](#stage-1--input-parsing)
- [第二階段 — 樣本正規化](#stage-2--sample-normalization)
- [第二階段 — 以 InChIKey 去除重複](#stage-2--deduplication-by-inchikey)
- [第二階段 — DAM（差異累積代謝物）](#stage-2--dam-differentially-accumulated-metabolites)
- [第三階段 — 富集分析（過度代表分析）](#stage-3--enrichment-over-representation-analysis)
  - [路徑模式](#pathway-mode)
  - [模組模式](#module-mode)
- [雙模式（正離子 + 負離子）輸入](#dual-mode-positive--negative-ionization-input)
- [快取與來源溯源](#caches-and-provenance)
- [儲存與載入工作階段設定（再現性）](#saving-and-loading-session-settings-reproducibility)
- [回報問題](#reporting-bugs)
- [主要參考文獻](#key-references)

---

<a id="stage-1--input-parsing"></a>
## 第一階段 — 輸入解析（Input parsing）

- **MS-DIAL `.txt`。** 前 4 列是 MS-DIAL 的中繼資料（`Class`、`File type`、`Injection
  order`、`Batch ID`）；第 5 列是欄位標題。當某一欄的 `File type` 值非空白、且不是
  `"NA"`、也不是字面上的列標籤 `"File type"` 時，該欄會被視為真正的樣本進樣——並保留
  在 `sample_cols` 中。這**包含** `Sample` 與 `Blank`（製程空白）；只排除 MS-DIAL 各
  組別的 `Average` / `Stdev` 彙總欄（標記為 `NA`）。被排除的樣本面積欄會顯示在底部面板
  的 **Data** 分頁（各槽位的輸入摘要），讓使用者隨時看到哪些欄被捨棄。
- **版本相容性。** 同時支援 MS-DIAL 4 與 MS-DIAL 5 的 Alignment 匯出檔。欄位是依名稱
  查找，因此 MS-DIAL 5 重新排序／改名的評分欄（它把 `Dot product` 拆成 `Simple` /
  `Weighted dot product`）也能以相同方式解析；metabolopan 只使用兩個版本共有的欄位。
- **缺失值。** 空白／僅含空白字元／`"null"`／`"NA"`／無法解析的強度儲存格會變成
  `f64::NAN`。明確寫成 `"0"` 的則維持 `0.0`。這與 `pandas.read_csv` 的語意一致，可避免
  下游統計把「缺測」與「真正的零」混為一談。
- **分組對應 `.csv`。** 標頭必須是 `sample,group`（嚴格兩欄）或 `sample,biosample,group`
  （三欄；`biosample` 欄為雙模式所必需——詳見下方 *雙模式輸入*）。其後的任何欄位都會被
  解析為選用的中繼資料。空白的 `group` 儲存格或重複的 `sample` 名稱會以明確的錯誤訊息
  拒絕。中繼資料欄會在載入時逐欄分類：若某欄的非空白儲存格全部都能解析為數字，就會出現
  在第二階段的「Metadata column」正規化單選按鈕中（例如 `dry_weight`、`dilution`、
  `total_protein`）；若某欄有任何非空白且非數字的儲存格（例如像 `CTR-01` 這樣的
  `biosample` 標籤），則會被靜默地排除在該單選按鈕之外，並在應用程式內的日誌窗格以一行
  WARN 說明被略過的是哪一欄、以及有多少儲存格無法解析。空白的中繼資料儲存格會解析為
  `None`。出現在 MS-DIAL `.txt` 中、卻不在 CSV 裡的樣本會被標記為 `Unassigned`（未指
  派）；CSV 中指名了 `.txt` 所缺樣本的列則會記錄為警告並忽略。**未指派樣本只在第一階段
  可見**——輸入摘要面板會以黃色的 `Unassigned (N samples)` 列顯示它們，讓你知道它們存
  在，但當你在第二階段設定畫面按下 `Start DAM` 時，它們就會從工作矩陣中被捨棄。正規化、
  去重複、DAM 統計或任何下游匯出都不會看到它們。若要讓某個樣本納入分析，請在中繼資料
  CSV 中為它加上真正的群組標籤；若要完全排除某個樣本（連第一階段都不顯示），請從
  MS-DIAL `.txt` 的 File type 列移除它那一欄（把該欄設為 `NA`）。第二階段下拉選單中的
  中繼資料欄順序會與 CSV 標頭順序一致（而非依字母排序），所以使用者看到的欄位順序就是
  當初書寫的順序。
- **第一階段 → 第二階段關卡。** `Continue to DAM` 按鈕會維持停用，直到：兩個檔案都成功
  解析、第 #1 槽位的離子化模式單選按鈕已設定、存在 ≥ 2 個不同的非 `Unassigned` 群組、
  且每個可指派群組都有 ≥ 2 個樣本（下游統計所必需）。此關卡是**與分析模式無關的**——
  分析模式（Analysis Mode）與 KEGG 物種／群組是稍後在第三階段設定畫面才設定，而非在此。

<a id="stage-2--sample-normalization"></a>
## 第二階段 — 樣本正規化（Sample normalization）

在進行任何逐特徵統計之前，使用者可選擇一種*樣本軸*（逐欄）正規化，以校正樣本之間的技術
性變異（進樣體積、稀釋倍數、乾重、總離子流）。每次 DAM 執行開始時，矩陣都會從最初解析得
到的 `intensity_raw` 重新正規化一次；`intensity_raw` 永遠不會被更動，因此切換方法是無損
的。預設為 `None`，會逐位元保留先前的行為。

除預設值外，另提供五種方法：

- **Sum（總和）。** 每個樣本的係數 = 該樣本所有非 NaN 強度的總和。輸出
  `x'[i, j] = x[i, j] / sum_j × median_j(sum_j)`。乘上各樣本總和的中位數可保留整體量級，
  讓 Welch / Student 路徑中選用的 `arcsinh` 步驟（由第二階段的 `Log transformation` 核取
  方塊控制；預設開啟）維持在有用的範圍內。
- **Median（中位數）。** 形式相同，改以各樣本的 NaN 感知中位數作為係數。
- **Metadata column（中繼資料欄）。** 使用者從中繼資料 CSV 解析出的選用數值欄中挑一欄
  （例如 `dry_weight`、`dilution`）。每個儲存格會除以該欄對應該樣本的值，再以所有樣本的
  中位數值重新縮放。資料不完整時的行為：
  - *缺值（空白儲存格）：* 該樣本會從分析中被**捨棄**——該樣本欄的每個儲存格都會標記為
    NaN，讓 DAM 的 NaN 感知機制將其排除在逐特徵統計之外。第二階段設定畫面會在使用者按下
    「Start DAM」之前，以一行黃色警告列出將被捨棄的樣本。
  - *非正值（零或負值）：* 會明確報錯並指出有問題的樣本與欄位。零／負的中繼資料是資料
    輸入問題、而非「缺值」，因此立即失敗（fail fast）才是正確做法。
  - *非數字儲存格：* 在 CSV 載入時即解析，並在抵達第二階段之前就報錯。
  - *群組前置檢查：* 在進行任何正規化工作之前，執行器會檢查：在捨棄沒有值的樣本之後，
    所選的分子組與分母組是否仍各自至少保有 2 個樣本。若否，錯誤橫幅會指出失敗的群組、
    欄位、剩餘數量、以及所需的最小值（`2`）。
- **Quantile（分位數）。** 強制讓每個樣本的分布對齊到一個共同的參考（各秩位在所有樣本
  間的平均）。本軟體遵循 Bolstad 與 Smyth 在 Bioconductor 支援討論串 #1569 中所達成共識
  的**原則**（2003，<https://support.bioconductor.org/p/1569/>）：在已排序位置
  `[k, k+t)` 出現並列（tie）的項目，會被指派為這 `t` 個秩位上參考值的**平均（MEAN）**
  ——`mean(reference[k..k+t])`。這是他們所達成的理論共識，而我們的程式碼如實地實作它。
  廣為部署的標準實作（`preprocessCore::normalize.quantiles`、`limma::normalizeQuantiles`）
  兩者皆採用較簡化的近似——平均秩查表搭配線性內插——僅在 `t == 2` 的並列、或參考在局部
  呈線性時才等於該原則，而在參考曲度較大時、對 `t ≥ 3` 的並列就會偏離該原則（這在以低於
  偵測極限值填補的代謝體學樣本底部很常見，例如下方的計算範例）。因此我們的輸出在這種情況
  下會與 preprocessCore / limma 不同。計算範例：參考為 `[1.5, 7.5, 52.5, 502.5, 55000]`，
  在已排序位置 1–3 有三項並列，在此會得到 mean(7.5, 52.5, 502.5) = **187.5**；
  preprocessCore / limma 則回傳 ref[2] = **52.5**。當所有樣本的非 NaN 數量相同時，兩派
  會產生完全相同的數字；差異純粹發生在參考不等長時。
- **Quantile — 各樣本非 NaN 數量不等。** 當樣本擁有不同數量的非 NaN 儲存格時（例如缺測
  情形不一致），參考會建立在大小為 `K = max(n_j)` 的共同分數秩格點上，並將每個樣本已排序
  的值線性內插到該格點上（與 limma 的 `(r − 1)/(n − 1) ∈ [0, 1]` 機制一致）。這可避免
  「較長的樣本主導高秩」這個錯誤——過去一個僅有 3 個非 NaN 的樣本，其最大值會被對應到
  參考的第 60 百分位（其 5 個位置中的 `reference[2]`），而非參考的第 100 百分位。當所有
  樣本的非 NaN 數量同為 `K` 時，每個分數秩都會落在整數格點索引上，內插路徑會退化為直接
  查表，輸出與本次變更之前我們所發布、僅支援等長的版本逐位元相同。NaN 儲存格維持 NaN。
- **PQN（Probabilistic Quotient Normalization，機率商數正規化）。** Dieterle 2006：先在
  內部做總和正規化；從所選群集（預設為 `All samples`，亦可選擇限定於某個指定群組）建立
  逐特徵的參考光譜；對每個樣本，計算其逐特徵商數相對於參考的中位數（略過參考為零、NaN、
  或樣本值為 NaN 的特徵）；再除以該係數並重新縮放。未指派樣本永遠不會抵達此階段（它們在
  第一階段 → 第二階段邊界就被捨棄了，因此無論是參考群集或逐樣本係數迴圈都不會看到它
  們）。若某個*已指派*的樣本仍產生退化的商數中位數（NaN 或 0），PQN 會中止並以
  `PqnDegenerateSamples` 列出有問題的名稱——請改用其他正規化方法，或從 MS-DIAL `.txt`
  的 File type 列移除該樣本。分派器的 INFO 日誌行會顯示一個 `reference_features_used=N`
  欄位（2026-05-29 新增），讓你看到該群集實際錨定了多少個特徵作為 PQN 參考（亦即
  `median(cohort) > 0` 者）相對於總特徵數——在不重跑流程的情況下，這對診斷 QC 稀疏度很
  有用。

**為何 Sum / Median / Metadata 是重新縮放到中位數係數（而非除到某個常數）。** 這三種方法
共用同一個機制。對每個樣本欄 *j*，它會計算一個純量係數 `f_j`——欄總和（Sum）、欄的 NaN
感知中位數（Median）、或樣本的正值中繼資料（Metadata）——然後把每個有限儲存格改寫為

    x'[i, j] = x[i, j] / f_j × M,    where   M = median_j(f_j)

`× M` 這一項是刻意的設計。總和正規化的教科書形式是單純的 `x / f_j`（或對 CPM 式計數乘上
`× 10^6`），這會強制把每個樣本拉到*每單位*尺度：對 Sum 而言，每欄總和都會變成 1（比例，
約 1e-5…1e-3）；對 Median 而言，每欄中位數都會變成 1。我們改為乘回 `M`，即**各樣本係數
的中位數**，使每欄的總和（或中位數）落在 `M`——*典型*樣本的原始量級——而非 1。樣本之
間的技術性負載（進樣體積、稀釋、乾重）仍被均衡；只有絕對強度尺度被保留下來。

- *為何重要——下游的 `arcsinh`。* 預設的 `Log transformation` 是 `arcsinh`，它只有在 *x*
  夠大時才表現得像對數（`arcsinh(x) ≈ ln(2x)`）；對接近 0 的 *x* 而言，它基本上是**線性
  的**（`arcsinh(x) ≈ x`）。把資料除成比例會把整個工作矩陣推進那個近線性區，使 arcsinh
  的變異數穩定化效果崩潰，並讓 t 檢定退化成在線性尺度上比較一堆極小的數。把數值維持在強度
  尺度（約 1e4–1e7）能讓 arcsinh 停留在它有用的類對數區——即「永遠不要接近 0」的目標。
  對**所有**儲存格乘上相同的常數 `M`，在 Brunner–Munzel 的中位數比值中會抵銷，但在
  `arcsinh` 之下**不會**（它是非線性的），因此這個重新縮放正是為了保護 Student / Welch +
  `arcsinh` 這條路徑——也就是目前的預設。
- *為何取係數的中位數（而非平均數）。* 中位數較穩健——單一個負載特別高的樣本無法把目標
  尺度往上拉——而且它讓*典型*樣本成為錨點：該樣本的 `f_j ≈ M`，於是 `f_j / M ≈ 1` 幾乎
  不會改變它，而偏離尺度的樣本則往它靠攏。
- *數字範例。* 三個樣本的欄總和為 `6, 15, 24`，得 `M = median(6, 15, 24) = 15`；做完
  `x / sum_j × 15` 後，每欄總和都變成 **15**（樣本 A 放大 ×2.5、C ×0.625、B 不變）——
  而非變成 1。若以各樣本中位數 `2, 20, 200` 做 Median 正規化，則 `M = 20`，每欄的中位數
  都會變成 **20**。Metadata 也相同，只是 `f_j` 改為所選欄的值，得到「在中位數乾重下的
  強度」。

與 MetaboAnalyst 對照：它的 Sum/Median 選項是除到固定常數（或比例），並搭配 `log10`
（後者還需要一個偽計數才能撐過零值）；metabolopan 的「先除再重新縮放到中位數」搭配的是
`arcsinh`，使正規化與廣義對數轉換在數值上保持相容。所選的 `M` 會在分派器 INFO 日誌中以
`scaling_to_median_factor=…` 回報。

**生命週期。** 正規化的選擇——以及其他每一個設定參數——會在整個工作階段的生命週期內、
跨越每一次導覽轉換而保留。退回上一階段絕不會丟失你的選擇；你只是回到上一個畫面，先前所有
的選擇都原封不動。（若你在第一階段重新挑選檔案，而先前的分子／分母組在新的中繼資料中已不
存在，第二階段會卡住關卡，直到你重新選擇有效的群組。）第三階段沒有獨立的正規化步驟——
第三階段富集分析看到的，就是那份（已正規化的）工作矩陣。

**啟動時的錯誤。** 正規化會在 DAM 的 tokio 任務生成之前同步執行，因此任何失敗（例如
`Sample 'A2' is missing a value in metadata column 'dry_weight'`）都會立即顯示在紅色橫幅
上。只有當工作矩陣為有限值且形狀正確時，DAM 任務才會啟動。

**值得知道的注意事項。**

- *Quantile* 假設各樣本的分布*理應*相同。對於同一基質、重複數充足的研究（例如細胞萃取
  物）這是合理的，但對於跨組織或跨生物體的比較——生物本質上在分布層級就有差異——則不
  成立。
- *PQN* 對大多數 NMR 式的稀釋變異很穩健。所選的參考群集很重要：當研究有一個乾淨的基線
  群組時，以它作為 PQN 參考，往往比 `All samples` 產生更清晰的生物訊號。**PQN 對樣本品質
  很嚴格**：若某樣本的逐特徵商數中位數為 `NaN`（沒有可用於對照參考的特徵）或 `0`（其非
  參考零特徵中有半數以上恰為 0——通常是稀疏／類空白的樣本），會以錯誤訊息列出有問題的
  樣本名稱。請從中繼資料 CSV 移除這些樣本，或改用較寬容的方法（None / Sum / Median /
  Metadata / Quantile）。在 2026-05-26 之前，這些樣本會在其他樣本被 PQN 縮放時靜默地維持
  未正規化，導致因尺度不一致而產生偏差的差異豐度判定。
- *Metadata* 的值必須為嚴格正值——除法與量級保留步驟都假設正值。零與負值會報錯，而非
  靜默地通過。
- *Sum/Median* 會完全保留樣本內的特徵比值；它們是同一種轉換的「縮放到量級」版本。兩者的
  差別在穩健性：Sum 對每個樣本中少數高強度離群值敏感；Median 則忽略它們。

<a id="stage-2--deduplication-by-inchikey"></a>
## 第二階段 — 以 InChIKey 去除重複（Deduplication by InChIKey）

MS-DIAL 經常為解析到同一個化合物的多個 Alignment ID 各自輸出一列。這有三種生物／儀器層面
的成因：

1. **加成物多重性。** 同一個中性分子在正離子模式下會以 `[M+H]+`、`[M+Na]+`、`[M+NH4]+`、…
   形式離子化（或在負離子模式下以 `[M-H]-`、`[M+Cl]-`、`[M+FA-H]-`、…）。每個加成物都會
   產生自己的 Alignment ID，但共用同一個 InChIKey。
2. **同位素峰。** MS-DIAL 會為 M0 單一同位素峰、以及 M+1 / M+2 天然豐度同位素峰各自輸出
   獨立的列（由 `Isotope tracking weight number`、或 `Adduct type` 中的 `[M+1]` / `[M+2]`
   後綴標示）。
3. **層析峰分裂。** 當峰偵測不夠理想時，單一個高斯沖提峰可能被切成兩個相鄰的 Alignment
   ID，它們共用每一項鑑定資訊，只在 `Fill %` / `S/N average` 上有差異。

把所有重複都餵進 DAM，會讓 FDR 的家族大小相對於真實化合物數膨脹 2–5 倍，侵蝕統計檢定力。
第三階段 ORA 的 `K`（抽取化合物數）也會膨脹，使路徑與模組的超幾何 p 值偏向那些恰好含有
「常見多加成物化合物」的條目。

**去重複以「預設啟用、可關閉」的切換開關呈現於第二階段設定畫面（預設開啟）。** 此階層判定
*純粹*是去重複作業，並非通用的品質過濾器——`inchikey = None` 的特徵會原封不動地通過，而
單一條目（一個 InChIKey 只對應一個 Alignment ID）即使鑑定品質不佳也會被保留。

### 階層判定表（Cascade decision table）

在每個相同 InChIKey 的群組內，存活的特徵由此階層中第一個能區分兩者的層級決定：

| 層級 | 欄位 | 判定規則 |
|-------|------------------------------------|-------------------------------------------------------------------------------------------------------|
| 1a    | `MS/MS matched`                    | `True` > `False` > 空白                                                                              |
| 1b    | `Total score`                      | 數值大者勝出（廠商計算的加權綜合分數，涵蓋所有光譜相似度指標，含 dot product） |
| 2     | 加成物類別                         | `Primary` > `NonPrimary` > `Dimer` > `Isotope`；在 `Primary` 之中，`[M+H]+` / `[M-H]-` > `[M+Na]+` / `[M+NH4]+` / `[M+K]+` / `[M+Cl]-` |
| 3a    | `Fill %`                           | 數值大者勝出（各樣本峰覆蓋率）                                                                |
| 3b    | `S/N average`                      | 數值大者勝出                                                                                           |
| 4     | `Alignment ID`                     | 字典序較小者勝出（決定性的最終判定）                                             |

加成物分類是決定性的且區分大小寫：`Isotope` 由 `Isotope tracking weight number > 0`、或
加成物字串中的 `[M+<n>]` 後綴偵測得出；`Dimer` 由開頭的倍率（`[2M+H]+`、`[3M-H]-`、…）
偵測得出；`Primary` 是封閉的允許清單 `{[M+H]+, [M+Na]+, [M+NH4]+, [M+K]+, [M-H]-,
[M+Cl]-}`；其餘一切（包含缺少加成物儲存格的情況）皆為 `NonPrimary`。

### 稽核 CSV（Audit CSV）

當 DAM 執行是在啟用去重複的情況下產生時，底部面板的 **Data** 分頁會在第二階段結果畫面上
顯示一個「Download dedup audit (CSV)」按鈕（在任何富集分析之後的畫面則不顯示）。CSV 格式：

```
# Deduplication audit — generated by metabolopan
# Total dropped: <N>; total kept: <M>; null-InChIKey passthrough: <K>
dropped_alignment_id,inchikey,winner_alignment_id,decided_at,loser_value,winner_value
```

`decided_at` 欄會告訴你每一次捨棄是由哪個階層層級所決定（`MsmsMatched` / `TotalScore` /
`AdductClass` / `FillPercent` / `SnAverage` / `Tiebreak`）；`loser_value` 與 `winner_value`
則承載決定性欄位在兩側各自的內容（若該側為 `None` 則留空）。在雙模式執行中，檔案會包含每
個模式各一份報告，以 `# Mode: POS` / `# Mode: NEG` 標頭行分隔。

### 關閉去重複（Opt-out）

在第二階段設定畫面取消勾選「Deduplicate features by InChIKey」即可停用。在未勾選的情況下，
DAM 執行與導入本功能之前的行為逐位元相同——每一個輸入列都會抵達前置過濾、FDR 的 `m`
等於前置過濾後的數量、且 DAM 結果上的 `dedup_report` 為 `None`。

<a id="stage-2--dam-differentially-accumulated-metabolites"></a>
## 第二階段 — DAM（差異累積代謝物，Differentially Accumulated Metabolites）

每個特徵都會在使用者所選的分子組與分母組之間被獨立檢定。本軟體提供三種統計方法；它們都遵循
相同的整體流程。

**1. 未知特徵過濾（預設開啟）。** `InChIKey` 為 `null` 的特徵（MS-DIAL 的「Unknown」鑑定；
約佔典型 Alignment Result 的 ~25%）會在任何統計工作之前被捨棄，這樣 FDR 校正的 `m` 就不會
納入那些終究無法進入第三階段 ORA 的代謝物。若使用者特別想對未鑑定特徵取得統計結果（例如標
記出供後續鑑定的候選），可在第二階段設定中取消勾選「Drop unknown features (no InChIKey)」
核取方塊。這是**相對於 Python 參考版本的一項偏離**——後者會把 Unknown 特徵留在 DAM 中、僅
在 PubChem CID 步驟才捨棄——這個對使用者可見的切換開關保留了以 Python 風格執行的選項。

**2. 逐特徵前置過濾。** 對每個剩下的特徵，會先捨棄合併 `numerator ∪ denominator` 欄中的
NaN 值，然後依序要求：(i) 分子組有 ≥ 2 個非 NaN 值、(ii) 分母組有 ≥ 2 個非 NaN 值、(iii)
合併後的 `nunique > 1`、以及 (iv) 合併後的 `IQR > 0`。未通過任一檢查的特徵會從結果中移除，
並計入 UI 中可見的 `skipped` 計數。檢查 (i) + (ii) 於 2026-05-29 新增——在此日期之前，某個
群組整組皆為 NaN 的特徵仍會通過前置過濾，並讓 NaN 經由檢定在結果中以一個 NS 槽位浮現；現在
它們會改在前置過濾層被略過，因此無法檢定的特徵不再於下游佔據 NS 槽位。

**3a. 方法：Student t 檢定（等變異數）** [參數方法，**預設**]。古典（同質變異數）形式。當
各組樣本數相近、且兩組離散程度大致相當時最適用——在這些假設下，它比 Welch 略具檢定力。
**新工作階段的預設值**：搭配 `Log transformation`（arcsinh）步驟（同樣預設開啟），它是本專案
的標準起點。若你懷疑變異數不等，請改用 Welch；若分布偏斜到連轉換都不足以應付，請改用
Brunner–Munzel。
- 與 Welch 共用的選用前置檢定轉換：當第二階段設定的 `Log transformation` 核取方塊被勾選時
  （預設開啟；`SessionSettings.log_transform = true`），會對每個非 NaN 儲存格套用
  `arcsinh(x)` 作為變異數穩定化步驟（asinh 能處理零／負值，而 log10 會把它們變成 NaN）。
  未勾選時，此步驟會被略過，原始工作矩陣的值會直接流入 t 檢定。較早寫死的 Pareto 縮放步驟
  已由 `add-log-transform-and-scaling`（封存於 2026-05-27）移除，因為實證已驗證逐特徵的線性
  重新縮放會在 t 統計量中抵銷——在 `log_transform=true` 之下，Welch / Student 的 p 值與變更
  前的流程逐位元相同。
- 合併變異數 `sp² = ((na − 1)·va + (nb − 1)·vb) / (na + nb − 2)`，固定自由度
  `df = na + nb − 2`，雙尾 p 值經由 Student-*t* CDF 求得。與 scipy 的
  `ttest_ind(equal_var=True)` 一致。
- **倍數變化的尺度取決於 `log_transform`。** 原因：在 `log_transform=true` 之下，*t* 統計量
  是在 arcsinh 轉換後的尺度上計算，但 arcsinh 對正值是凹函數，因此由 Jensen 不等式可知，兩
  個重尾群組的*原始*平均比值，可能在**正負號**上與 *t* 檢定實際評估的 arcsinh 平均差不一致。
  若把原始平均比值與 arcsinh 尺度的 *p* 值並列回報，會靜默地誤判由離群值驅動的特徵（例如
  `num=[0.1]×9 + [100]` vs `den=[5]×10`，得原始 FC ≈ 2.02 ⇒「Up」，但 Welch *t* ≈ −3.25、
  *p* ≈ 0.01 ⇒「Down」）。依 `arcsinh-scale-fc-on-log-transform` 變更（2026-05-29），參數
  方法分支現在會讓尺度一致：
    - `log_transform=false`（原始尺度）：`FC = mean(numerator) / mean(denominator)`，
      `log2(FC) = log2(FC)`。`FcBasis::Mean`。與 2026-05-29 之前的流程逐位元相同。
    - `log_transform=true`（arcsinh 尺度）：`log2(FC) = (mean(arcsinh(num)) −
      mean(arcsinh(den))) / ln(2)`，且 `FC = 2^log2(FC)`。`FcBasis::ArcsinhMean`。在相同資料
      上，`log2(FC)` 的正負號**保證**與 *t* 統計量的正負號一致。對大的 *x*，arcsinh(x) ≈
      ln(2x)，因此 `log2(FC)` 會漸近於 `log2(GM(num) / GM(den))`——即 limma / DESeq2 的古典
      對數倍數變化。對小的 *x*（接近 0），arcsinh(x) ≈ x，因此 `log2(FC)` 會退化為縮放後的
      算術平均差，而非真正的比值。這是變異數穩定化已被記載的結果；等價的對數 FC 詮釋只在大
      *x* 漸近區（arcsinh 與 ln 對齊處）成立。
  CSV 匯出會經由 `fc_basis` 欄（`mean` / `median` / `arcsinh-mean`）標示目前作用的基準，讓
  下游使用者無需重跑流程即可辨識某個數字位於哪種尺度上。

**3b. 方法：Welch t 檢定（異變異數）** [替代的參數方法]。與 Student 屬於同一參數族，但不
假設變異數相等。當各組離散程度明顯不同時、或當你不確定而想用較安全的預設時，請用它。
- 與 Student 相同的選用前置檢定轉換（僅 `arcsinh`，由第二階段 `Log transformation` 核取方塊
  控制；預設開啟）。Pareto 縮放已由 `add-log-transform-and-scaling`（2026-05-27）對兩條參數
  路徑一併移除。
- Welch 的 t 統計量是用（可選地經 arcsinh 轉換的）值、以 NaN 感知的平均數與變異數計算，搭配
  Welch–Satterthwaite 自由度，再經由 Student-*t* CDF 轉換為雙尾 p 值。與 scipy 的
  `ttest_ind(equal_var=False)` 一致。
- **倍數變化的尺度與檢定尺度一致**——規則同上方 Student。在 `log_transform=true` 下，`FC`
  位於 arcsinh 尺度（`FcBasis::ArcsinhMean`），故其正負號永遠與 Welch *t* 的正負號一致。在
  `log_transform=false` 下，`FC` 是古典的原始平均比值（`FcBasis::Mean`）。

**Welch / Student 的邊界情況——某一組變異數為零。** 當某一組的每個重複樣本都有相同的值
（例如該特徵在某個條件的每個樣本中都低於偵測極限、而被填補為一個常數）時，
Welch–Satterthwaite 自由度會塌縮為*另一*組的 `n − 1`。對 `n = 2` 而言這給出 `df = 1`，使
*t* 分布變得極寬、p 值非常保守——即使兩組在視覺上明顯分離也是如此。這是標準的數學行為（與
R 的 `t.test(var.equal=FALSE)` 與 SciPy 的 `ttest_ind(equal_var=False)` 完全一致），但對代謝
體學而言，受影響的特徵往往對應到你可能想保留的「一個條件有、另一個條件無」的真實訊號。
`run_dam` 每次執行會發出單一行 INFO 日誌，回報觸發此路徑的特徵數（在你的工作階段日誌
`<binary>/data/logs/session_*.log` 中尋找 `zero_variance_features=N`）；當 N > 0 時，可考慮
改用 Brunner–Munzel 方法重跑，它以秩為基礎，對此邊界情況的處理方式不同。自 2026-05-29 起，
此診斷計數器採用相對容差——變異數低於 `(max(|mean|, 1))² × 1e−20` 即被標記——因此那些其群組
在浮點數雜訊範圍內為常數的特徵（例如在高強度尺度上對逐位元相同的正規化前輸入做運算，其
`var ≈ ε² × c²` 雖非零、但 df 病態仍以相同方式發作）也會計入。t 檢定函式內部各方法自有的
`var == 0.0` 防護則維持不變——只放寬了這個診斷計數器。

**3c. 方法：Brunner–Munzel + Cliff's δ** [無母數方法]。當各組的強度分布偏斜或不等、且變異
數穩定化轉換仍不足時適用。代謝體學資料常難以用高斯假設充分描述（高度偏斜的對數分布、頻繁的
有無模式、批次假影），在這些情況下，Brunner–Munzel + Cliff's δ 能在此工作流程所見的各種離散
情形下提供更誠實的 p 值。當預設的 Student t 檢定（即使在 `arcsinh` 之後）並不適配時——例如
高度偏斜或以有無為主的特徵、或為了對應先前已發表的無母數分析——請經由第二階段設定的單選按鈕
選用它。
- Brunner–Munzel 統計量是以 `numerator ∪ denominator` 上的中位秩計算，結合類
  Welch–Satterthwaite 自由度，再經由 Student-*t* 分布轉換為雙尾 p 值。行為與 SciPy 的
  `brunnermunzel(distribution='t')` 與 R 的 `lawstat::brunner.munzel.test` 一致——`sqrt` 內
  的 W 分母為 `nx·Sx + ny·Sy`（2026-05-26 之前的實作使用 `(nx+ny)·(Sx/nx + Sy/ny)`，在等 n
  時會把 `|W|` 膨脹 `sqrt(N/2)` 倍；修正前在每組 n=5 時，BM p 值系統性地約偏顯著 1.58 倍。
  先前以 BM 產生的 CSV / 火山圖匯出應重新產生）。
- Cliff's δ 效應量：`(gt − lt) / (n · m)`，其中 `gt` 與 `lt` 分別為「嚴格大於」與「嚴格小
  於」的配對計數。範圍 −1 .. +1；此處採用慣用的「中等效應」門檻 |δ| ≥ 0.33。
- 倍數變化使用各組**中位數**：`FC = median(numerator) / median(denominator)`，且
  `log2(FC) = log2(FC)`。中位數對離群值穩健，與以秩為基礎的檢定哲學一致。

**4. 多重檢定校正。** 每次第二階段執行都會對逐特徵 p 值套用使用者所選的 FDR 校正，無論這些
p 值是由哪種統計方法產生。第二階段設定畫面提供一個含兩個選項的單選按鈕：

- **Benjamini–Hochberg（BH）**——預設值。與 R `p.adjust(method='BH')` 及 MetaboAnalyst 的
  慣例一致，使第二階段結果可直接與那些工具已發表的數字比較。BH 假設各檢定之間獨立、或具正
  迴歸相依。
- **Benjamini–Yekutieli（BY）**——選擇性啟用。將 BH 的 q 值乘上精確的調和級數因子
  `c(m) = Σ_{i=1}^{m} 1/i`（對大 m 約為 ln(m) + γ，故在 m = 5,000 時 BY 大約比 BH 保守 9
  倍）。BY 在任意正相依下都能控制 FDR，因此當許多特徵在生物學上彼此相關時（例如共享同一路徑
  成員的代謝物），它是較安全的選擇。

NaN 的 p 值在任一方法下校正後皆仍為 NaN。所選方法會回報於火山圖的註解列上（例如
`FDR(BH)<0.05`），並以開頭的 `# FDR: BH` / `# FDR: BY` 註解行寫入每一份 DAM CSV 匯出，使螢幕
截圖與下載檔皆能自我說明。參考文獻：Benjamini & Hochberg (1995)；Benjamini & Yekutieli
(2001)。

**5. 趨勢分類**（隨使用者調整門檻值即時重算——從不存入結果）。預設門檻值：`FC = 2.0`（等同
|log2(FC)| ≥ 1.0）、`FDR = 0.05`、`|δ| ≥ 0.33`（僅 BM）。
- Student / Welch（皆為參數方法，無效應量）：當 `FDR < threshold` 且
  `log2(FC) > log2(fc_threshold)` 時為 `Up`；當 `FDR < threshold` 且
  `log2(FC) < −log2(fc_threshold)` 時為 `Down`。δ 門檻對參數檢定會被忽略。
- Brunner–Munzel：上述參數規則**且** `|δ| ≥ delta_threshold`。`δ = None` 的特徵（BM 因某一
  組少於 2 個非 NaN 值而無法計算效應量）會被分類為 `NotSignificant`。

**6. 火山圖。** X 軸 = `log2(FC)`，Y 軸 = `−log10(p_adjusted)`。**X 軸代表什麼，取決於目前作
用的方法與 `Log transformation` 切換開關**——對 `log_transform=false` 的 Welch / Student 是
平均比值、對 `log_transform=true` 的 Welch / Student 是 arcsinh 平均差（以 log2 為單位）、對
Brunner–Munzel 是中位數比值。詳見上方第 3 節；目前作用的基準記錄於每個 `DamFeature` 的
`fc_basis`（`mean` / `arcsinh-mean` / `median`）。三種顏色對應趨勢分類（紅 / 藍 / 灰，α ≈
0.5）。門檻線為黑色虛線：在 `−log10(FDR)` 的水平線，以及在 `±log2(FC)` 的垂直線。`log2(FC)`
為 `±∞` 的特徵（某一組的平均或中位數恰為 0）會被停靠在 X 軸邊緣 `±(xabs_max + 0.5)`，並加上
小幅抖動以保持可見。Y 軸的對稱飽和處理：BH/BY q 值下溢為恰好 `0.0` 的特徵（極大的 `|t|` / 極
小的原始 p，常見於分離良好的群組）會被停靠在 Y 軸頂端（`y_max`）**正下方**，每點向下抖動至多
`0.08`（以 `−log10(q)` 為單位；與 X 軸 ±0.04 抖動慣例的尺度相符），以免多個飽和特徵堆疊在同
一像素上（2026-05-29 新增——在此日期之前，所有 q=0 的特徵都對應到恰好 `y_max` 而在視覺上塌
縮）。底層的 `neg_log10_p_adjusted` 值仍為 `f64::INFINITY`，並如此記錄於 CSV 匯出——只有畫面
上的位置被抖動。Y 軸在其他情況下僅為了顯示而裁切於 `finite_max + 1`；底層數值仍保留於 CSV
匯出。`NaN` 的 `neg_log10_p_adjusted` 保留給真正「p 無法計算」的情況（BM 完全分層的群組；參
數檢定在 NaN 捨棄後 `n < 2`）——這些點會從圖中捨棄，但仍列於 CSV。X 軸標籤下方有一條註解
列，摘要說明方法、目前作用的 FC 基準（`FC: mean` / `FC: median` / `FC: arcsinh-mean`，
2026-05-29 新增，使 X 軸的詮釋無需查閱 CSV 即可明確）、目前門檻值、以及 ±∞ 計數——例如
`Method: Brunner-Munzel | FC: median | FDR(BH)<0.05, FC≥2.0, |δ|≥0.33 | −∞: 12  +∞: 8`。

**BM 點的大小編碼 Cliff's δ 的量值。** 在 Brunner–Munzel 的繪圖上，每個散布點的半徑由該特徵
的 `|Cliff's δ|` 對應而來：`|δ|=0` 給出最小但仍可見的點、`|δ|=1` 給出約 1.3 倍預設半徑的點，
中間量值則在兩個錨點之間線性縮放。右側圖例會在既有的趨勢計數下方長出第二個 `|δ| size` 區
塊，含三個位於 `|δ|=0/0.5/1.0`、以中性灰色呈現的參考點——把散布點與這些參考點做大小比對即可
從圖上讀出量值。Welch / Student 的繪圖在整張圖上維持一致的點半徑，且**不**繪製 `|δ| size`
圖例區段（那些檢定不會產生可供編碼的 Cliff's δ）。`|δ|` 未定義的 BM 特徵（某一組 `n < 2` 個
非 NaN 值）會退回預設半徑，並仍以對應的趨勢顏色繪製。

**第二階段值得知道的注意事項。**
- BM 以中位數為基礎的 FC，意味著小 n 研究（例如每組 3 個樣本）比 Welch 以平均為基礎的 FC 更
  容易產生 `±∞` 的 log2(FC)，因為三個樣本中只要有一個零就會把該組中位數拉到零。註解列會顯示
  ±∞ 計數，故此情況絕不會是靜默的。
- 每組 n = 2 的參數 t 檢定（Student 或 Welch）只有約 1–2 個自由度、並不可靠；第一階段關卡的
  「每組 ≥ 2 個樣本」要求讓你不致低於下限，但檢定力充足的參數檢定希望每組 ≥ 5。當等變異數假
  設成立時，等樣本數的 Student 是三者中最靈敏的；當變異數明顯不同時，Welch 是穩健的備案。
- 趨勢分類取決於目前作用的門檻值。CSV 匯出器寫入的是匯出當下所計算的趨勢，與火山圖在該時刻
  所顯示的完全相同。

<a id="stage-3--enrichment-over-representation-analysis"></a>
## 第三階段 — 富集分析（過度代表分析，over-representation analysis）

第三階段接收第二階段的 DAM 結果，並提問：*「在我這份差異豐度化合物清單中，哪些 KEGG 條目
被過度代表？」*——這裡的「條目」指的是**一條 KEGG 路徑**（路徑模式）或**一個 KEGG 模組**
（模組模式）。兩種模式在超幾何檢定、使用者所選的 FDR（BH 或 BY）、以及可測代謝體母體上共用
完全相同的機制；它們只在 ORA 所操作的化合物集合目錄上有所不同。

第三階段擁有自己獨立的 FDR 校正單選按鈕，與第二階段的選擇無關——預設為 BH 以利跨工具再現
性，但**對路徑／模組 ORA 而言，建議選擇 BY**，因為條目之間本質上會共用化合物（許多化合物
出現在多條路徑中），這違反了 BH 的獨立性假設。第三階段點圖的色標標題與註解列都會標明目前
作用的方法（例如 `-log10(FDR (BY))` / `FDR: BY`），而富集分析 CSV 開頭的 `# FDR:` 註解行也
會記錄此選擇供下游解析。

### 富集分析設定畫面（Enrichment Analysis setup screen）

第三階段設定畫面是使用者進行下列選擇之處：

- **分析模式（Analysis Mode）**（Pathway / Module），以單選切換鈕選取。兩種模式的選擇**以及**
  其各自擷取的 KEGG 快取，會在整個工作階段的生命週期內並存——在模式間切換是即時的，絕不會
  重新擷取你已經拉取過的資料。
- **KEGG 範圍。** 路徑模式顯示一個可搜尋的物種選擇器，內含預先載入的 KEGG 生物清單；模組模式
  則顯示下方 *模組模式* 所述的 Level + Group 選擇器。選定一個物種（或 Group）會在此畫面就地
  觸發對應的 KEGG 擷取——一個附帶說明文字的小型進度列會串流顯示逐路徑（或逐模組 + ETA）的
  進度，無須離開設定畫面。
- **要納入的方向（Direction）**（`Both` / `Up only` / `Down only`）。
- **某路徑／模組中偵測到的最少化合物數**（即「最小條目大小」過濾器；標籤會隨模式調整——路徑
  模式為「… in a pathway」、模組模式為「… in a module」；預設 `1`，範圍 `[1, 20]`）。會在
  建立 FDR 家族**之前**，捨棄其**母體限定**化合物數（`m_p = |entry.compounds ∩ universe|`）
  低於此門檻的路徑／模組。`m_p`（此處）與超幾何的 `m` 參數都使用交集的**集合**基數：在某個
  KEGG 條目的 COMPOUND 區塊中被列出超過一次的化合物只**計一次**，而非按出現次數計——按原始
  出現次數計會膨脹 `m`（以及期望命中數）並壓低該條目的富集比（一個於 2026-05-26 修正的重複
  計數錯誤）。預設 `1` 很寬鬆——只有 `m_p = 0` 的條目會被捨棄（它們本來就會短路為
  `p = 1.0`），因此每條至少含一個可測化合物的路徑都會被檢定。把它調高到 **`3` 以對應
  MetaboAnalyst 的 `minPathSize`**，可排除那些雖小但在數學上無法檢定的 `m_p ∈ {1, 2}` 條目
  ——`m_p = 1` 的條目無論命中多強都永遠無法產生 BH 顯著的結果，因此捨棄它們可減少多重檢定的
  懲罰（代價是檢定的路徑較少）。與 *最小命中數* 正交：這個旋鈕是會縮小 `m` 的 **FDR 前條目
  過濾器**；*最小命中數* 則是不改變 p 值的 **FDR 後顯示過濾器**。
- **富集分析 FDR 門檻**（預設 `0.05`）。
- **FDR 校正**（預設 BH；ORA 建議 BY——詳見上文）。
- **最小命中數**（FDR 後顯示過濾器；預設 `1`）。
- **`Run Enrichment`** 按鈕（擷取進行中時停用；停用狀態的滑鼠懸停提示會說明是哪個擷取正在
  阻擋按鈕）。控制點圖顯示上限的 Top N 輸入欄位位於下一個畫面（Enrichment Result），讓你可以
  在看過資料後再迭代調整，無須回到設定畫面。

<a id="pathway-mode"></a>
### 路徑模式（Pathway mode）

流程為：

1. **身分解析（PubChem PUG REST）。** 對每個通過第二階段前置過濾的特徵（**不只是** DAM 顯著
   者），經由對 `compound/inchikey/property/InChIKey/CSV` 的 POST，把其 `InChIKey` 解析為一
   個或多個 PubChem CID。每批最多 200 個 InChIKey。
2. **KEGG 化合物轉換（KEGG REST）。** 對每個唯一的 CID，經由
   `/conv/compound/pubchem:CID1+CID2+...` 解析為一個 KEGG 化合物（`cpd:Cxxxxx`）。每批最多
   10 個 CID，由 KEGG 用戶端節流（每次請求間隔 334 ms，約 3 req/s，符合 KEGG 文件所載的軟性
   上限）。對應到 `glycan:` 或 `dr:` 的列會被濾除——只保留 `cpd:` 目標。HTTP 403 會被視為速率
   限制訊號，並以 5 秒退避重試至多 5 次。
3. **多重對應規則。** 一個特徵就是一個化學實體。若 PubChem 為某個 InChIKey 回傳多個 CID（立
   體／區域／鹽類的歧義）、而它們全部解析到同一個 KEGG cpd，則該特徵只把該 cpd **計一次**到
   DAM 化合物集合 `K` 與母體 `N` 中。若它們解析到確實不同的 cpd（對正常的 PubChem 紀錄而言
   極為罕見），則每個 cpd 各自計入。`InChIKey` 沒有 PubChem CID、或其 CID 全部無法對應到 KEGG
   cpd 的特徵，會從 `K` 與 `N` 中捨棄，並顯示於底部面板 **Data** 分頁的對應漏斗中
   （`<N> InChIKeys → <N> PubChem CIDs → <N> KEGG cpds`）。
4. **母體定義（N）。** 母體是所有「通過第二階段前置過濾**且**成功經 PubChem 與 KEGG conv 對應
   成功」的已鑑定特徵之唯一 cpd ID 聯集——即此平台上的*可測代謝體*。這比 MetaboAnalyst 的
   「此物種中所有 KEGG 化合物」母體更保守，後者會用你的儀器偵測不到的化合物來膨脹 `N`。我們
   刻意採用「僅可測」母體，使 p 值更能反映你的資料原本能說明的範圍。
5. **FDR 前的條目大小過濾。** 在任何超幾何工作之前，每條路徑的 `M_p` 會與使用者可調的
   `min_entry_size`（預設 `1`，範圍 `[1, 20]`）比較。`M_p < min_entry_size` 的條目會從本次
   執行中**完全捨棄**——它們不對 FDR 家族貢獻任何 p 值、不出現在 CSV 中、也不出現在點圖上。
   被捨棄的數量會顯示於底部面板 **Data** 分頁的一行保留率資訊
   `Tested: <surviving> / <total> (≥ N compounds in universe)`
   （模組模式為 `Tested: <surviving> (≥ N compounds in universe)`）。預設 `1` 讓前置過濾保持
   寬鬆——只有 `M_p = 0` 的條目會被捨棄（它們本來就會短路為 `p = 1.0`），因此每條至少含一個
   可測化合物的路徑都會被檢定。把它調高到 **`3` 以對應 MetaboAnalyst 的 `minPathSize`**，並
   額外排除在典型 `K`/`N` 值下數學上無法檢定的 `M_p ∈ {1, 2}` 條目——例如 `M_p = 1` 的條目
   最多只能產生 `k_p = 1`，得到原始 `p ≈ K/N`，這很少低於 `α = 0.05`，更少低於 BH 臨界值
   `0.05/m`。取捨是對稱的：較低的 `min_entry_size` 會檢定更多路徑，但也擴大多重檢定家族 `m`。
6. **逐路徑超幾何檢定。** 對每條通過條目大小過濾而存活的路徑 `p`，以
   `M_p = |unique(pathway.compounds) ∩ universe|`（該路徑落在可測母體內之唯一 cpd ID 的集合
   基數——單一 COMPOUND 區塊內的重複 cpd ID **不會**膨脹 `M_p`）與
   `k_p = |K ∩ pathway.compounds|`：
   - `p_value = 1 - HypergeometricCDF(k_p - 1; N, M_p, K)`（看到「至少」 `k_p` 個命中的上尾
     機率）
   - 若 `k_p, M_p, K, N` 中任一為零，實作會短路為 `p_value = 1.0`（避免未定義的 CDF 參數）。
   - 實作還以 `debug_assert!` 強制要求 `K ⊆ N`（任何讓 K 漏出 N 之外化合物的上游迴歸，都會在
     開發／測試中被大聲地捕捉；發行版建置會發出每次執行的 `ERROR` 日誌，摘要說明任何
     Hypergeometric 定義域錯誤，使「所有條目皆不顯著卻無任何診斷」這種失敗模式無法靜默出貨）。
7. **使用者所選的 FDR 校正**，經由第三階段設定畫面的獨立單選按鈕（預設 Benjamini–Hochberg；
   Benjamini–Yekutieli 一鍵可及；**None** 作為第三個選項，僅供探索性執行，詳見下文）。此單選
   按鈕刻意與第二階段的選擇相互獨立（依封存的 `add-bh-fdr` design.md D3）：兩個階段有不同的
   相依性樣態，許多使用者會合理地想要第二階段 BH（火山圖的跨工具再現性）+ 第三階段 BY（對共享
   化合物條目採保守 ORA）。
   對路徑／模組 ORA，即使 MetaboAnalyst 與多數生物工具預設 BH，我們仍**建議 BY**：路徑大量
   共享化合物（糖解 ↔ TCA 共享 G6P、丙酮酸等），因此 BH 所依據的獨立性假設被違反。BY 在相依
   下較保守；在相同資料上，預期得到比 MetaboAnalyst 一致更高（較不顯著）的校正後 p 值。
   **None** 會完全略過多重檢定校正——結果表與 CSV 中的 `fdr` 欄會原封不動地承載原始 p 值。
   請**僅用於探索性排序**，絕不可用於已發表的顯著性宣稱；在典型的 KEGG 路徑目錄（約檢定 300
   條路徑）上，純粹出於機率，你會預期在 `p < 0.05` 處有約 15 個偽陽性。第二階段 DAM 的單選
   按鈕**不**提供 None——對約 13k 個特徵的原始 p 會淹沒結果集；手工製作、攜帶
   `dam_fdr_method=NoCorrection` 的快照會被防禦性地強制改回 BH，並伴隨一個 `tracing::warn!`
   事件。
   **色階。** 每個標記的填色編碼 `-log10(FDR)`（在 None 之下為原始 `-log10(p)`），採用
   ColorBrewer **YlOrRd** 9 階漸層——所顯示中最不顯著的條目（FDR 位於顯示門檻處）為最淡的
   黃色，往最顯著者加深為深紅；點與色標圖例共用單一個 `-log10` 跨距，因此相同顏色在兩者間代表
   相同顯著性。
   目前作用的方法會記錄於點圖的色標標題（`-log10(FDR (BH))` / `-log10(FDR (BY))` / None 時為
   `-log10(p-value)`——外層包裝被去除，因為軸上的值「就是」原始 p、而非 q），並記錄於匯出之
   富集分析 CSV 開頭的 `# FDR: BH` / `# FDR: BY` / `# FDR: None` 行。CSV 還攜帶第二行註解
   `# MinEntrySize: N`，記錄該次執行所用的 FDR 前過濾門檻，使檔案能自我說明。點圖本身在 X 軸
   下方還附帶一個四行的純文字註解區塊，讓審閱者僅憑圖即可重建 FDR 家族：
   `Background universe = <N> compounds measured and mapped to KEGG` /
   `Compounds of interest = <K> differentially abundant (increased | decreased | both directions)` /
   `Pathways tested = <m>[ of <total>  ·  <dropped> skipped (< <min_entry> compounds each)][; ≥ <min_hit> hits required]` /
   `Significance: FDR-adjusted, Benjamini–Hochberg (BH)`（或 `… Benjamini–Yekutieli (BY)`、
   或 `raw p-value (no FDR correction)`）。`N` / `K` / `m` 這些符號刻意以文字拼寫而非縮寫；
   被檢定數 `<m>` 是抵達 BH/BY 的條目數，也是每個原始 p 值被乘上的除數。
   `m` 這個分母等於**通過 FDR 前 `min_entry_size` 過濾**（步驟 5）而存活的路徑數——亦即
   `m = entries.len() − entries_dropped_by_min_entry_size`。在預設 `min_entry_size = 1` 下，
   `M_p = 0` 的前置過濾條目也會被捨棄。協調器層級的 Group 過濾（模組模式）在更早的一層套用；
   到 FDR 執行時，`m` 已反映了這兩道過濾。
8. **顯示過濾（FDR 後）。** 一個使用者控制的 `min_hit_count`（預設 1）會把命中數較少的路徑從
   點圖與 CSV 中隱藏。這是*顯示*過濾器——`m` 早已在所有存活條目上計算完畢，因此無論此設定為
   何，FDR 值都是誠實的。與步驟 5 的 `min_entry_size` 不同：那個是會縮小 `m` 的 **FDR 前條目
   過濾器**；這個是不改變 p 值的 **FDR 後顯示過濾器**。
9. **點圖的「選取」與「排序」——兩種不同的依據。** 點圖以**刻意不同的準則**來決定*哪些*條目
   被繪出、以及*如何*把它們在 Y 軸上堆疊：
   - **選取（哪些條目出現）依統計顯著性。** 在通過 `fdr < threshold` 與 `min_hit_count` 過濾
     （步驟 7–8）的條目中，圖會保留 **FDR 最低的 Top N**（`top_n`，預設 20，可在結果畫面調
     整）。因此所顯示的條目永遠是*最顯著*的那些——它們**絕不**依富集倍數選取。
   - **垂直順序（Y 軸）依效應量。** 被保留的條目接著依**富集倍數（觀測／期望）由大到小**排列，
     使富集倍數最大的條目位於**最上方**，整張圖沿 X 軸（X 軸本身即為富集倍數）讀起來像一道
     「大者在上」的階梯。同分時先依 FDR（較顯著者在前）、再依條目 ID 決定。這符合
     clusterProfiler「以 X 軸度量排序 Y 軸」的慣例。
   實務上的後果：當顯著條目多於 `top_n` 時，被略去的是**最不顯著**（FDR 最高）的那些——*而非*
   富集倍數最小的那些。顯著性把關「能否納入」；效應量只負責排列「已納入者」。匯出的 CSV 與此
   無關：它列出每一個存活條目，依 FDR 由小到大排序，並附上完整（未截斷）的名稱。
10. **點圖畫布高度。** 匯出圖的高度會自動配合實際顯示的列數——`clamp(min(top_n, displayed) ×
    0.3 + 1.0, 2.0, 40)` 英吋——並在你每次 Draw / Re-draw 時**重新計算**。因此若某次執行在你
    最初的 FDR 門檻下不顯著，而你在結果畫面放寬門檻後重繪，畫布會長大以容納新顯現的列，而非
    把它們塞進一張短圖（那會截斷 Y 軸標籤）。編輯 **Height (in)** 欄位會把它變成手動覆寫，並
    一直維持到下一次富集分析執行／重跑重置自動配合為止。

    **文字大小與條目數無關。** 標籤、軸標題、色標、以及 Hits 圖例都隨圖的**寬度**縮放（固定
    的 `Width (in) × DPI`），*而非*自動配合的高度——因此兩個條目的結果，其文字渲染大小與二十
    個條目的結果完全相同。（在此修正之前，稀疏的圖會以短畫布為依據決定字型大小，結果文字又小
    又難讀。）高度的 `2.0` 英吋下限——從 `1.5` 提高而來——存在的目的是讓全尺寸圖例在那些稀疏
    結果上永遠能容於畫布之內。

<a id="module-mode"></a>
### 模組模式（Module mode）

模組模式執行與路徑模式完全相同的 PubChem → KEGG conv → 超幾何 → 使用者所選 FDR 流程，但
**(a)** 條目目錄是 KEGG 模組的集合、而非逐物種的路徑，且 **(b)** 使用者挑選的是一個**生物
群組（organism Group）**、而非單一物種。當某模組的 KEGG `COMPLETE` 區塊含有所選 Group 中至少
`min_group_overlap` 個生物時，該模組才會被納入分析；這就是逐物種框架對應到全域模組目錄的方式。

1. **生物群組選擇。** 當分析模式切換鈕設為 Module 時，第三階段 **Enrichment Analysis setup**
   畫面會顯示一個 Level 單選按鈕（1 / 2 / 3）與一個 Group 下拉選單。Level 對應到 KEGG
   `/list/organism` 的譜系欄（Level 1 為 `Eukaryotes`、Level 2 為 `Animals` / `Bacteria` 等、
   Level 3 為 `Mammals` / `Insects` 等）。KEGG 目前有 11,744 個生物，全部恰好有 4 個譜系層級；
   我們公開前三層。挑選一個 Group 會具現化 `org_codes`：屬於該 Group 的 KEGG 生物代碼集合
   （`hsa`、`ath`、…）。（在封存於 2026-05-25 的 `reorder-gui-and-move-mode-to-stage3` 變更
   之前，此選擇器位於第一階段；現在它位於第三階段設定，讓你只在看過 DAM 結果之後才確定模式。）

2. **模組 → Group 過濾。** 每個模組的 `/get/<module-id>` 回應攜帶一個 `COMPLETE` 區塊，列出
   該模組被完整組裝的生物。當以下條件成立時，模組會被保留供 ORA：
   ```
   |module.complete_orgs ∩ group_orgs|  ≥  min_group_overlap
   ```
   預設 `min_group_overlap = 1` 很寬鬆（∃-重疊：Group 中任一個生物即足夠）。較高的值會收緊
   過濾——例如 `min_group_overlap = 5` 要求 Group 的生物中至少有 5 個完整組裝了該模組。目前
   作用的門檻會顯示於第三階段結果標頭，使你發表的任何數字僅憑「標頭 + 快取快照」即可再現。

3. **母體與 K——同路徑模式。** PubChem 與 KEGG-conv 階段與模式無關。`N` 仍是可測代謝體（對應
   到某 KEGG cpd 的 DAM 特徵）；`K` 仍是符合目前作用方向過濾（Up / Down / Both）之 DAM 特徵的
   cpd 集合。模組模式*不會*用「所有模組化合物」或「所有 KEGG 化合物」替換 `N`。

4. **逐模組超幾何檢定。** 與路徑模式相同：對每個被保留的模組 `m`，
   `M_m = |module.compounds ∩ universe|`、`k_m = |K ∩ module.compounds|`，且
   `p_value = 1 - HypergeometricCDF(k_m - 1; N, M_m, K)`，採用相同的零輸入短路。

5. **使用者所選的 FDR 校正**——選項與預設皆同路徑模式（預設 BH；對共享化合物條目建議 BY）。
   `m` 分母等於（經 Group 過濾後）**被保留的模組數**，而非 KEGG 目錄中全部約 573 個模組。這是
   正確的虛無：ORA 問的是「在那些*可能*適用於此生物群組的模組中，哪些被過度代表？」把分類學上
   不相關的模組納入 `m`，會在不貢獻生物訊號的情況下把 FDR 往上扭曲。

6. **空 COMPOUND 模組計數器。** 有些 KEGG 模組（標誌型／僅含反應的模組）根本沒有 `COMPOUND`
   區塊。它們的條目會以 `compounds = []` 通過 ORA 並短路為 `p_value = 1.0`。底部面板 **Data**
   分頁會以一行 `With compound list: <kept>  (−<empty> empty)` 顯示這些，使靜默捨棄絕不會侵蝕
   信任。（對等的路徑模式回報已列入規劃。）

**模組模式值得知道的注意事項。**
- **首次執行的成本。** 從 KEGG 冷擷取目前列出的全部約 573 個模組，在 334 ms 的請求間節流
  （3 req/s）下約需 6–8 分鐘。模組 ID 範圍為 `M00001`–`M01063`，但 KEGG 讓此範圍保持稀疏
  ——已退役的 ID 不會重新使用，因此實際數量低於上界。第三階段設定畫面會顯示一個行內進度列，
  其 ETA 在最初 5 個模組完成後，由逐模組實際耗時的滾動平均推導而得。後續執行會使用快取，
  `Run Enrichment` 按鈕在數秒內即可啟用。
- **Group 1 只有兩個選項**（Prokaryotes / Eukaryotes），在生物學上非常粗略。它的存在是為了
  完整性——例如「任一原核生物」的比較研究——但多數分析會受益於 Level 2（6 個候選）或 Level 3
  （數十個候選）以取得更精細的範圍界定。
- **`min_group_overlap` 是一個研究旋鈕。** 預設 `1`（寬鬆的 ∃-重疊）適合探索性工作。對論文而
  言，可考慮測試較高的門檻以確保穩健性——一個 Group 中數百個生物（例如「Animals」）裡只有一
  個擁有的模組，對該分析框架而言在生物學上是邊緣的，即使它通過了預設過濾。
- **模組 CSV 欄名與路徑模式 CSV 一致。** 兩種模式都匯出相同的標頭：
  `EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,FDR,HitKeggIDs`。在模組模式
  中，`EntryID` 欄攜帶 `M00001` 式的模組 ID；在路徑模式中則攜帶
  `<species_code><pathway_number>` 形式的 ID（例如 `gmx00010`）。

### 開始新一輪分析（Starting a new analysis round）

當你完成一次富集分析、想分析另一份資料集——或從頭重跑整個流程——時，第三階段 **Enrichment
Result** 畫面會在 `[Download enrichment results CSV]` 下方獨立一行提供一個 **Start a new
analysis** 按鈕。按下它會開啟一個確認對話框，警告目前的 DAM / 富集分析結果、以及任何尚未下載
的圖或 CSV 都將遺失。按下 **Start over** 後，應用程式會把每個參數重設為其預設、清除已載入的
MS-DIAL `.txt` / 中繼資料 `.csv` 以及記憶體中的 KEGG 資料，並把你帶回第一階段——*但不會*
重跑啟動時的生物清單載入。（磁碟上的 KEGG 快取仍保留，因此之後重新擷取相同物種或模組會是快速
的快取命中。）

這刻意有別於階段步進器的 **Input** 步驟，後者會導覽回第一階段，同時*保留*每一項設定、已載入
的檔案、以及已擷取的快取，讓你能在**同一份**資料集上持續迭代。用步進器來微調並重跑目前分析；
用 **Start a new analysis** 來捨棄一切、從頭開始。若你之後可能還想用到目前的設定，請在開始
重來之前，透過 Data 分頁的 **[Save settings…]** 按鈕儲存它。

<a id="dual-mode-positive--negative-ionization-input"></a>
## 雙模式（正離子 + 負離子）輸入（Dual-mode input）

代謝體學實驗常把同一批生物樣本同時跑正離子化與負離子化兩種模式，每個研究因此產生兩個
MS-DIAL `.txt` 匯出檔。本應用程式支援一次載入這兩個檔案，並在一條保守的聯集規則下合併它們的
富集訊號。

### 何時使用雙模式

只要你對同一批生物樣本同時擁有 POS 與 NEG 兩個 `.txt`、並想要一份能反映「任一離子化所提供
證據」的單一富集結果，就使用雙模式。單模式（一個 `.txt`）仍為預設。（自 2026-06-01 起，單模式
的 K 也套用「僅衝突即嚴格（conflict-only-strict）」規則——詳見下方*第三階段 — 雙模式 N 與 K
的運算*——因此它與先前版本逐位元相同，**唯獨**對於「某化合物同時被一個 Up 特徵與一個 Down
特徵抵達」的資料集而言，該化合物現在會被排除。）

### 準備輸入

1. **兩個 `.txt` 檔案。** 每個離子模式各一。`Adduct type` 欄同時驅動：槽位 1 模式單選按鈕的
   自動填入（見下方*第一階段 UI*），以及當使用者手動覆寫為相反極性時的一個建議性「不一致」
   提示（以 `+` 結尾的加成物推論為 Positive，以 `-` 結尾的推論為 Negative）。
2. **一個 3 欄的中繼資料 CSV**，標頭為 `sample,biosample,group`。每一列把一個逐模式的樣本名稱
   （例如 `CTR_positive_01`、`CTR_negative_01`）對應到其**生物樣本標籤**（兩個模式皆為同一個
   `CTR-01`）與群組。biosample 欄讓工具能辨識兩個名稱不同的樣本其實是同一個生物重複。

以 2 欄 `sample,group` CSV 進行的雙模式執行，會在第一階段被一個明確的錯誤擋下——請新增
`biosample` 欄、或移除第二個 `.txt` 以繼續。

> **單模式並不需要 `biosample` 欄。** 它只在載入第二個 `.txt` 時才為必需。只有一個 `.txt`
> 時，單純的 `sample,group` 形式就足夠（若存在 `biosample` 欄，會以名稱辨識並排除於第二階段
> 中繼資料正規化單選按鈕之外——它不會被當作數值中繼資料欄提供）。

### 不平衡或缺少某模式的樣本

下方的計算範例是完全平衡的（每個生物樣本都在兩個模式中跑過），但真實研究有時只在單一極性下
採集某個生物樣本。`biosample` 欄是配對兩個模式的依據，因此第一階段會在允許 `Continue to DAM`
**之前**強制執行三項雙模式完整性檢查。每一項都會顯示一個明確的錯誤：

1. **每個群組在*每個*模式中都需要 ≥ 2 個樣本。** 若某群組在 POS 有足夠的重複、但在 NEG 掉到
   2 以下（例如因為數個生物樣本缺少 NEG 採集），第一階段會擋下並顯示
   `Group 'X' has N sample(s) in POS but M in NEG — both modes need ≥ 2.`。只要每個群組仍各自
   達到「每模式 2 個」，少數缺某模式的樣本是可以容忍的；把關的只有「逐群組、逐模式」的計數。
2. **生物樣本在同一模式內必須唯一。** 兩列把同一個生物樣本標籤對應到同一個模式，會觸發
   `Biosample 'B' appears in N POS rows — must be unique per mode.`。
3. **生物樣本在跨模式間必須維持同一群組。** 若 `CTR-01` 在 POS 是 `control`、在 NEG 卻是
   `treatment`，第一階段會擋下並顯示 `Biosample 'B' is in group 'X' in POS but 'Y' in NEG.`。

**通過關卡的缺某模式樣本所造成的影響。** 兩個模式各自在自己的樣本欄上獨立執行 DAM——某個在
NEG 缺席的生物樣本，單純就不會在 NEG 執行中被迭代，因此該模式對其群組的重複較少、檢定力也
相應較低；但這不會使該次執行失效。在第三階段，聯集是在**化合物**層級建立（依下方「僅衝突即
嚴格」規則）、而非樣本層級，因此某生物樣本缺少一個模式，只會讓該模式對受影響的化合物貢獻
`Absent`——整合後的 K 不受影響。

**建議。** 為了得到最乾淨的雙模式結果，請在兩種極性下都採集每個生物樣本。若某些樣本確實只有
單一極性，要嘛只在它們存在的那個模式中保留它們（只要每個群組仍各模式 ≥ 2）、要嘛捨棄不平衡的
那一側。出現在某 `.txt`、卻不在中繼資料 CSV 中的樣本會被標記為 `Unassigned`，並在第一階段 →
第二階段邊界自動捨棄（見上方分組對應的說明），這是排除不想要欄位的另一種方式。

### 第一階段 UI

槽位 #1（永遠可見）與槽位 #2（由 `+ Add second ionization mode` 按鈕揭示）各有一個檔案選擇
器、一個模式單選按鈕（Positive / Negative）、以及一份逐槽位摘要。在 `auto-infer-stage1-ion-mode`
（封存於 2026-05-26）之後，槽位 1 的模式單選按鈕會在每次新檔載入與重新挑選時，自 `infer_polarity(&table)`
自動填入：`≥ 95%` 為正離子後綴的 Adduct 欄會將其設為 Positive、`≥ 95%` 為負離子後綴設為
Negative、模稜兩可的混合則讓單選按鈕保持未設定（既有的灰色「Could not auto-detect…」提示仍
適用）。當槽位 1 的模式已設定時，槽位 2 的單選按鈕會在三種觸發下自動填入**相反**值：(1) 槽位
2 經由 `+ Add second ionization mode` 按鈕揭示、(2) 槽位 2 的 `.txt` 被載入、(3) 槽位 1 的模式
改變（手動點擊或重新挑選重新推論）——情況 (3) 在新的槽位 1 值與槽位 2 已顯示的值衝突時，也會
翻轉槽位 2。使用者仍可手動點擊任何單選按鈕加以覆寫。槽位 2 的單選按鈕仍會停用已被槽位 #1 選走
的選項（以提示說明原因）。加成物不一致提示（「黃色：Adduct 欄說是 X，但你選了 Y」）仍會在與
自動偵測相牴觸的手動覆寫時出現；兩種提示都不會擋下 `Continue to DAM`。第一階段畫面**不**帶有
分析模式切換鈕或任何 KEGG 選擇器——那些位於 DAM 之後的第三階段設定畫面。

### 第二階段（共用設定，逐模式 DAM）

第二階段使用單一設定畫面——一種正規化方法、一組比較（分子組 vs 分母組）、一種 DAM 方法、一種
FDR 方法——並套用於**兩個**模式。在協調器內部，兩個 tokio 工作者會逐模式平行執行 `run_dam`；
執行中畫面會顯示兩條堆疊的進度列。若任一模式失敗，錯誤訊息會指出是哪個模式（`Positive: ...`
或 `Negative: ...`）。

火山圖畫面會在繪圖區上方渲染一條 `POS | NEG` 分頁列。每個分頁快取自己的材質；改變任一門檻滑桿
會使兩者都失效。PNG 匯出使用逐模式的預設檔名（`volcano-pos.png` / `volcano-neg.png`）。DAM 的
CSV 匯出會發出開頭的 `# Mode: dual (POS+NEG)` 註解行，並在每一列前面加上一個 `Mode` 欄，列的
順序為 POS 在前、NEG 在後。

### 第三階段 — 雙模式 N 與 K 的運算（Stage 3 — dual-mode N and K math）

第三階段在「僅衝突即嚴格」的聯集規則下，從**兩個**模式的 DAM 特徵建立母體 N 與前景 K（此為
保守的選擇：方向相反的訊號會排除該化合物）。

**PubChem 與 KEGG `/conv` 呼叫對聯集後的 InChIKey 集合只執行一次**，因此雙模式下的網路成本
不會加倍。

**N（母體）** = 經由 PubChem → KEGG conv 鏈、可從任一模式中任一特徵抵達的每個 cpd 之聯集。

**逐模式趨勢彙總。** 對每個 cpd `c`，分別從各模式蒐集逐特徵趨勢並彙總：
- `Up`       — 此模式中有任一貢獻特徵標記為 Up、且無 Down
- `Down`     — 對稱
- `NS`       — 只有不顯著特徵
- `Conflict` — 同一模式中同時有 Up 與 Down 特徵（同 InChIKey、不同趨勢的邊界情況）
- `Absent`   — 該 cpd 完全無法從此模式抵達

**「僅衝突即嚴格」規則下的 K（前景）。** 對目前作用方向 `Up`：當至少一個模式說 Up、**且**沒有
任何模式說 Down、**且**沒有任何模式處於 Conflict 時，該 cpd 才進入 K。`Down` 為對稱。`Both`
要求至少一個 Up 或 Down 訊號、且無 Conflict、且非（一個模式 Up 而另一個模式 Down）。

**單模式套用相同的衝突規則（自 2026-06-01 起）。** 單模式執行是此規則的退化單一模式情況：在
單一模式內被一個 Up 特徵與一個 Down 特徵同時抵達的化合物——兩個不同的 InChIKey 對應到**同
一個** KEGG 化合物、一個 Up + 一個 Down——會彙總為 `Conflict` 並從 K 中**排除**，這與雙模式
是相同的保守選擇。（在此之前，單模式會把這類有歧義的化合物留在 K 中。）被衝突排除的數量會
出現在第三階段 INFO 日誌中。對於沒有這類模式內衝突的任何資料集，單模式 K 維持不變。

底部面板 **Data** 分頁會把雙模式的分割顯示為母體／前景來源溯源漏斗的一部分：

```
Universe — all tested features (measurable metabolome)
  … → N KEGG cpds  (POS-only: a; NEG-only: b; in both: c)
Foreground — significant features (active direction)
  … → K KEGG cpds  (sig POS-only: d; sig NEG-only: e; agree both: f; excluded by conflict: g)
```

當某一個模式貢獻了全部的 K cpd 時，會出現一行黃色的
`K source: POS only (NEG had 0 sig features in the active direction)`。富集分析 CSV 會發出開頭
的 `# Mode: dual (POS+NEG)` 註解行；逐列的 CSV 形狀不變（ORA 的運算與模式無關）。

### 計算範例

使用 `data/double-mode/` 的測試固定檔（3 個 control + 3 個 treatment 生物樣本 × 2 個模式 =
12 個樣本欄 + 6 個生物樣本）：

1. 第一階段：把 `POS_*.txt` 載入槽位 #1（Mode: Positive）、把 `NEG_*.txt` 載入槽位 #2
   （Mode: Negative）、並載入 `metadata.csv`。按下 `Continue to DAM`。
2. 第二階段：選 `treatment` vs `control`，正規化與 FDR 維持預設。執行中畫面會顯示兩條進度列；
   視特徵數而定，每個模式約需 6–60 秒。
3. 第二階段門檻：在 POS 與 NEG 分頁間切換以檢視各自的火山圖；下載一個分頁式 PNG 或統一的
   CSV。按下 `Continue to Enrichment`。
4. 第三階段設定：挑選一個 KEGG 物種（路徑模式）或 Level + 生物群組（模組模式）；行內進度列
   會串流顯示 KEGG 擷取。完成後，按下 `Run Enrichment`。
5. 第三階段結果：結果面板會顯示分解區塊；被衝突排除的 cpd ID 會以 INFO 出現在日誌中。若你想要
   點圖上多一些／少一些列，就地調整 Top N，然後按下 `Re-draw dot plot`。點圖會保留 Top N 個
   *最顯著*的條目、並依*富集倍數*堆疊（大者在上）；畫布高度會在每次重繪時重新配合所顯示的列
   （見上方「點圖的選取與排序」）。

<a id="caches-and-provenance"></a>
## 快取與來源溯源（Caches and provenance）

第三階段的快取（`pubchem.json`、`cid_to_cpd.json`、`modules.json`）儲存**逐條目**的
`fetched_at: DateTime<Utc>` 時間戳，與第一階段物種快取（檔案層級時間戳）不同。逐條目的粒度是
刻意的：這些快取會在數週或數月、跨多個工作階段中逐步增長，而檔案層級時間戳要嘛謊報年齡、要嘛
強迫頻繁的全量重新整理。第三階段結果畫面會把這呈現為一個時間跨距
（`PubChem: 2026-03-01 -> 2026-05-22`）；模組模式還會額外顯示模組快取在該次執行所用之**被保留**
模組上的時間跨距，而非整個快取。

快取鎖：
- **PubChem `.inchikey.lock` + KEGG `.cid_to_cpd.lock`**——短暫存在，僅在快取寫入期間持有。
  30 秒等待、100 ms 輪詢。（兩個檔案都以點號開頭／為隱藏檔。）
- **KEGG `.modules.lock`**——長時間執行的建議鎖，在整個約 6–8 分鐘的模組擷取期間持有。鎖檔嵌入
  持有者的 PID 與一個心跳 `last_seen_at` 時間戳（至多每 30 秒改寫一次）。並行的應用程式實例看到
  活著的鎖時，會等待至多 30 分鐘（5 秒輪詢）直到它清除。若心跳超過 90 秒未更新，該鎖會被視為
  孤兒（持有者已當機）並被覆寫。這可避免兩個應用程式實例競相跑過模組擷取迴圈、並一起觸發 KEGG
  的 403 速率限制。
- **啟動時清理。** 每次應用程式啟動時，快取目錄中所有的 `.lock` 檔案都會被無條件移除，因此當機
  絕不會永久阻擋未來的寫入。

快取新鮮度——**無過期門檻**。所有 KEGG 快取都不會過期：無論年齡多大，被快取的條目永遠會被回
傳，且應用程式從不自行靜默重新擷取。（早期版本中的 7 天路徑 / 30 天模組門檻已於 2026-05 移除。）
取而代之的是，底部面板 **Data** 分頁的 `Cache data` 區塊（位於 Enrichment Analysis +
Enrichment Result 畫面）會中性地呈現擷取時間，把重新整理的決定權留給你（這些擷取時間行與
Refresh 按鈕已於 2026-05 從畫面主體移到 Data 分頁）：
- 逐物種路徑快取：顯示 `Cached <ts> (N days ago)`（設定畫面）或 `KEGG pathways (<code>): <ts>`
  （結果畫面）；經由 `Refresh KEGG pathway cache` 按鈕重新擷取。
- 模組條目：顯示一個 `KEGG modules fetched date: <oldest> -> <newest>` 跨距；暖擷取的決定取決
  於快取鍵的成員資格（此模組是否已被快取？），而非年齡。經由 `Refresh KEGG module cache` 按鈕
  重新擷取。
- 在 Enrichment **Result** 畫面，目錄重新整理按鈕（模組 / 路徑）會導覽回 Setup 畫面以在那裡執行
  重新擷取（其進度列位於該處）；PubChem / KEGG-conv 的重新整理則經由確認對話框就地執行。
- 生物清單（`organisms.json`）：在啟動時載入一次（快取優先：無論年齡多大，磁碟上的副本永遠勝
  出），且沒有應用程式內的重新整理按鈕。若要強制重新擷取一份新的 `/list/organism`，請從快取
  目錄刪除 `organisms.json` 並重新啟動。（`Refresh KEGG pathway cache` 按鈕只會重新擷取所選
  物種的「路徑→化合物」對應，而非整份名冊。）

<a id="saving-and-loading-session-settings-reproducibility"></a>
## 儲存與載入工作階段設定（再現性）（Saving and loading session settings）

Data 分頁中的兩個按鈕——**[Save settings…]** 與 **[Load settings…]**——讓你把每一個第一至
第三階段的參數快照到一個 JSON 檔，並於日後重新套用。其用意在於再現性：若你（或合作者）以相同
快照與相同輸入重跑，分析會逐位元相同。

### 檔案內容

一份美化排版的 JSON，包含：

- `schema_version`（目前為 `3`——`log_transform: bool` 先經 `add-log-transform-and-scaling`
  把它從 `1 → 2`，接著 `min_entry_size` 再經 `add-min-entry-size-filter` 把它從 `2 → 3`，兩者
  皆封存於 2026-05-27）、`app_version`、`saved_at`（UTC）、以及一個初始為 `""` 的 `user_note`
  欄位——你可以用任何文字編輯器打開檔案並填寫它。
- `input_files`——對你在儲存當下所載入的每個 MS-DIAL `.txt` 與中繼資料 `.csv`：該檔的基本檔名
  + 其 SHA-256。**只有雜湊值——你的原始資料絕不會被包含。** 這讓未來的 Load 能偵測你的輸入是否
  已偏離當初製作快照時的版本。
- `settings`——從第一階段到第三階段的每一個參數（分析模式、物種／生物群組、比較群組、DAM 方
  法、正規化、FDR 方法、門檻值、匯出尺寸、富集分析方向 / FDR / top-N）。

### 各按鈕何時可用

- **Save settings…** 在啟動啟始畫面之後的每個畫面上都啟用，無論是否已載入輸入。從空白的第一
  階段儲存，可把你偏好的預設值擷取為下次的預設組態。
- **Load settings…** **只在第一階段**啟用。在其他階段，按鈕呈灰色；懸停時顯示「Loading
  settings is only available on the Stage 1 input screen.」。這是刻意的——在分析中途套用快照會
  讓畫面上的結果與新參數不同步，因此此工作流程要求你從輸入重跑。

### 載入流程

1. 在第一階段點擊 **[Load settings…]**。作業系統的檔案選擇器會開啟。
2. 挑選一個已儲存的 `.json`。一個確認對話框會向你顯示其內容：
   - 儲存當下的時間戳（以你的當地時間）、快照的 app 版本、使用者備註（若有）。
   - 一行設定摘要（分析模式、DAM 方法 + FDR、正規化、富集分析方向 + FDR + top-N）。
   - **雜湊不符**——若你目前載入的任何輸入檔與快照的 SHA-256 不同，會在此列出。若你繼續，設定
     仍會套用，但你會被警告輸入已偏移。
   - **欄位重設**——若快照指名了某個分子／分母組、某個中繼資料欄、或某個 PQN 參考群組，而它在
     你目前載入的中繼資料中不存在，這些欄位會被列出並在套用時重設為 `None`。你需要在第二階段
     設定重新挑選它們。（此區段只在你於 Load 當下已載入中繼資料時才出現；若你在上傳中繼資料之
     前就 Load，安全網改由第二階段設定關卡承擔——見下一段。）
3. 點擊 **Apply settings** 以覆寫你目前的設定，或點擊 **Cancel** 以捨棄。

### 若我在上傳中繼資料之前就載入設定，會怎樣？

快照的 `numerator` / `denominator` 會被原封不動地寫入設定（Load 當下不做驗證，因為沒有中繼
資料可供比對）。當你之後上傳中繼資料並前進到第二階段設定時，關卡會檢查群組成員資格：若被保留
的值未出現在新中繼資料的群組中，「Start DAM」按鈕會呈灰色，並附帶一個行內警告
（`⚠ Numerator/denominator group not present in the loaded metadata.`）以及相同文字的懸停提
示。從 ComboBox 下拉選單重新挑選一個有效群組，警告即會清除。

### 手動編輯 JSON

此檔為純 UTF-8 JSON、美化排版。你可以：

- 在 `user_note` 欄位加上備註。
- 微調單一門檻值，而無須從應用程式重新儲存。
- 移除 `input_files` 區塊以分享「僅設定」的快照（Load 會處理空的 `input_files` 陣列——雜湊
  檢查會被略過）。

把 `schema_version` 手動編輯成 `3` 以外的數字、或破壞 JSON 語法，會在 Load 時顯示一個明確的
錯誤提示（例如 *"This settings file uses schema version 1; this app expects version 3."* 或
*"Settings file is not valid JSON (line 7 column 15) …"*）。在 `add-min-entry-size-filter`
（2026-05-27）之前儲存的快照攜帶 `schema_version == 1` 或 `== 2`，兩者都會被拒絕——請從你目前
的設定重新儲存，以產生一份 v3 快照。

<a id="reporting-bugs"></a>
## 回報問題（Reporting bugs）

若有什麼看起來不對勁——非預期的錯誤、某個卡住的階段、與預期不符的結果——取得協助最簡單的方式
是點擊日誌窗格（視窗底部、**Clear** 按鈕旁）中的 **[Download bug report…]**。一個確認對話框會
列出產生的 zip 將包含哪些檔案，接著一個存檔對話框讓你選擇要放在哪裡。

該 zip 恰好包含八個檔案：

- `README.txt`——說明此套件包及其隱私界線。
- `version.txt`——app 建置資訊（套件版本、git SHA、rustc、target）。
- `RUST_LOG.txt`——僅 `RUST_LOG` 指示詞的值，於單一行上。
- `KEGG_CACHE_DIR.txt`——僅 `KEGG_CACHE_DIR` 環境變數值（或 `<unset>`）。這兩者是逐變數的檔案
  （檔名 = 變數名），讓任何人都不會把此套件包誤認為完整的環境傾印——只有這兩個指名的變數會被
  納入。
- `logs.txt`——本工作階段中每一個 INFO / WARN / ERROR 事件（HTTP 與其他低階相依套件的雜訊會被
  濾除，以保持檔案易讀）。
- `app_state.txt`——你當時所在的階段與你目前的設定（分析模式、物種／群組、比較群組、FDR 方法、
  門檻值等）。
- `input_summary.txt`——你所載入的 MS-DIAL 檔案與中繼資料 CSV 的路徑與計數（僅路徑——無儲存格
  值）。
- `cache_summary.txt`——KEGG / PubChem 快取檔案的大小與新鮮度時間戳（無快取內容）。

**隱私：**

- 此套件包絕不包含你的原始 MS-DIAL `.txt` 輸入、你的中繼資料 CSV、或任何先前的 CSV/PNG 匯出。
- 套件包內的絕對路徑會把你的家目錄替換成 `~`（例如 `/Users/alice/Projects/study/POS.txt` 變成
  `~/Projects/study/POS.txt`），使套件包在公開分享時（GitHub issue、電子郵件）不會洩漏你的帳號／
  使用者名稱。
- 只有 `RUST_LOG` 與 `KEGG_CACHE_DIR` 環境變數會被呈現——絕不會是完整的處理程序環境。

你可以放心地把此 zip 附加到 GitHub issue 或以電子郵件寄出，無須擔心洩漏你的實驗資料或機器身分。

逐工作階段的日誌檔也會保存在磁碟上的 `<binary-directory>/data/logs/` 下，保存 7 天，之後於啟動
時自動刪除。若你想擷取先前某次執行的日誌，請在重新開啟應用程式之前先到該目錄查看。

<a id="key-references"></a>
## 主要參考文獻（Key references）

- **Brunner–Munzel 檢定。** Brunner, E. & Munzel, U. (2000). *The nonparametric
  Behrens-Fisher problem: Asymptotic theory and a small-sample approximation.* Biometrical
  Journal 42 (1): 17–25.
- **Cliff's δ。** Cliff, N. (1993). *Dominance statistics: Ordinal analyses to answer
  ordinal questions.* Psychological Bulletin 114 (3): 494–509.
- **Welch t 檢定。** Welch, B. L. (1947). *The generalization of "Student's" problem when
  several different population variances are involved.* Biometrika 34 (1/2): 28–35.
- **BY FDR。** Benjamini, Y. & Yekutieli, D. (2001). *The control of the false discovery
  rate in multiple testing under dependency.* Annals of Statistics 29 (4): 1165–1188.
- **BH FDR。** Benjamini, Y. & Hochberg, Y. (1995). *Controlling the false discovery rate:
  A practical and powerful approach to multiple testing.* Journal of the Royal Statistical
  Society B 57 (1): 289–300.
- **超幾何 ORA。** Khatri, P., Sirota, M. & Butte, A. J. (2012). *Ten years of
  pathway analysis: Current approaches and outstanding challenges.* PLoS Computational
  Biology 8 (2): e1002375.
- **KEGG。** Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. & Ishiguro-Watanabe, M.
  (2023). *KEGG for taxonomy-based analysis of pathways and genomes.* Nucleic Acids
  Research 51 (D1): D587–D592.
- **KEGG MODULE。** Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M. & Tanabe, M.
  (2014). *Data, information, knowledge and principle: back to metabolism in KEGG.*
  Nucleic Acids Research 42 (D1): D199–D205.
- **分位數正規化。** Bolstad, B. M., Irizarry, R. A., Åstrand, M. & Speed, T. P.
  (2003). *A comparison of normalization methods for high density oligonucleotide array
  data based on variance and bias.* Bioinformatics 19 (2): 185–193.
- **機率商數正規化（PQN）。** Dieterle, F., Ross, A., Schlotterbeck,
  G. & Senn, H. (2006). *Probabilistic quotient normalization as robust method to account
  for dilution of complex biological mixtures. Application in 1H NMR metabonomics.*
  Analytical Chemistry 78 (13): 4281–4290.
