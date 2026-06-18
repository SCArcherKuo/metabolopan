# 使用手冊

軟體版本：metabolopan v1.2.3

更新日期：2026-06-19

本手冊記錄本軟體在數值上的運作方式——演算法、預設門檻值、以及與常見替代做法的差異——讓你能在論文或報告中為它產生的每一個數字辯護。
在發表任何依賴本軟體的結果之前，請先閱讀一次本手冊。

本管線的運作分為三個階段。

<a id="how-to-read-this-manual"></a>
## 如何閱讀本手冊（How to read this manual）

本手冊採用**雙軌**書寫，因此無論你只是想跑一次分析，或是需要為論文中的每一個數字辯護，它都行得通。

- **平實語言導讀**以 2–4 句話為每個章節開頭：說明該步驟做什麼、你如何點擊操作、你會挑選什麼、又會得到什麼回饋。只要快速瀏覽這些導讀，你就能跑完整條管線。
- **更深入的技術區塊**（公式、邊角情況，以及任何標示 `> **給開發者：**` 的內容）位於其下方。初次閱讀時你可以跳過它們，等到審稿人問「為什麼是這個數字？」時再回頭查看。

**格式圖例**（全文一致套用）：

- **粗體** = 你點擊或設定的 UI 控制項——按鈕、單選按鈕、核取方塊或欄位（例如 **Start DAM**、**DAM method** 單選按鈕、**Log transformation** 核取方塊、**Continue to DAM**）。
- `等寬` = 在磁碟上或檔案中看到的字面文字——檔名、CSV 標頭、設定鍵、日誌行、公式與值（例如 `0.05`）。
- `> **給開發者：**` 提示框收錄非寫程式者永遠不需要的內部實作細節。可以放心跳過。

你還會看到其他提示框：`> **注意：**`（澄清，或本軟體與標準工具的差異）、`> **⚠ 警告：**`（資料輸入錯誤與快速失敗的條件），以及 `> **範例：**`（計算範例與直覺）。

**三條閱讀路徑：**

- **（a）只想跑一次分析。** 閱讀 [第一階段 — 輸入解析](#stage-1--input-parsing)、[第二階段 — 正規化、去除重複與 DAM](#stage-2--normalization-deduplication--dam)、[第三階段 — 富集分析（過度代表分析）](#stage-3--enrichment-over-representation-analysis) 的導讀，以及 [計算範例](#worked-example)。
- **（b）為論文中的數字辯護。** 閱讀方法子章節（[DAM](#differentially-accumulated-metabolites-dam) 檢定方法 3a–3c）、[多重檢定校正（FDR）](#4-multiple-testing-correction)、[缺失值 vs 真正的零](#missing-values-nan-vs-true-zeros-00)，以及 [主要參考文獻](#key-references)。
- **（c）再現性／撰寫指令稿。** 閱讀 [儲存與載入工作階段設定（再現性）](#saving-and-loading-session-settings-reproducibility) 與 [快取與來源](#caches-and-provenance)。

<a id="pipeline-at-a-glance"></a>
## 管線一覽（Pipeline at a glance）

**本軟體適合誰：** 濕實驗代謝體學研究者——不需要寫任何程式就能執行。

本管線的運作分為三個階段。每個階段接收一個輸入、執行一項關鍵操作，再把結果交給下一個階段：

| 階段 | 輸入 | 關鍵操作 | 輸出 |
|-------|-------|---------------|--------|
| **第一階段 — 輸入** | MS-DIAL `.txt` + 群組對應 `.csv` | 解析並對齊，讓*缺失*儲存格與*真正的零*保持區分 | 一張解析後的表格 |
| **第二階段 — 正規化 → 去除重複 → 檢定** | 解析後的表格 | 選用的樣本正規化，接著以 InChIKey 去除重複，再進行逐特徵統計檢定 | **DAM 特徵 + 一張火山圖** |
| **第三階段 — 富集分析** | DAM 特徵 | InChIKey → PubChem CID → KEGG 化合物，接著做超幾何過度代表分析 | **一張點圖** |

> **注意：** 單模式跑一個 MS-DIAL `.txt` + 一個群組 `.csv`；雙模式跑兩個 MS-DIAL `.txt` 檔（正離子 + 負離子）+ 一個含 `biosample` 欄的群組 `.csv`。你所選的模式在第一階段設定，並貫穿後續每一個階段。

<a id="key-terms--symbols"></a>
## 關鍵術語與符號（Key terms & symbols）

**符號**

| Symbol | 平實語言意義 | 階段 |
|--------|------------------------|-------|
| N | 母體大小——分析所能取用的 KEGG 化合物唯一總數 | 3 |
| K | 前景大小——被抽取（DAM）的唯一 KEGG 化合物 | 3 |
| m | FDR 校正所套用的特徵數量（檢定家族大小） | 2 |
| m_p | 逐路徑 ORA 的檢定家族大小 | 3 |
| k_p | 落在某個給定路徑中的被抽取化合物 | 3 |
| M_m | 某個給定模組目錄中的化合物 | 3 |
| k_m | 落在某個給定模組中的被抽取化合物 | 3 |
| M（中位數因子） | 各樣本正規化因子的中位數——每一欄被重新縮放所朝向的目標量級 | 2 |
| f_j | 樣本欄 *j* 的逐樣本正規化因子（欄總和、中位數或中繼資料值） | 2 |
| δ | Cliff's δ 效應量——兩群組分離的程度，從 `−1` 到 `+1` | 2 |
| FDR | 偽發現率——每個特徵經校正後的顯著性值（q 值） | 2 |
| FC | 分子組與分母組之間的倍數變化 | 2 |
| log2(FC) | 倍數變化的以 2 為底對數——火山圖的 X 軸 | 2 |
| q-value | 某個特徵經 FDR 校正後的 p 值 | 2 |
| c(m) | BY 乘在 BH q 值上的調和因子 `∑_{i=1}^{m} 1/i` | 2 |
| df | t / Brunner–Munzel 分布的自由度 | 2 |

**術語**

| Term | 一句話意義 |
|------|----------------------|
| DAM | 差異累積代謝物（Differentially Accumulated Metabolites）——在所選兩群組之間有顯著差異的特徵。 |
| ORA | 過度代表分析（Over-representation analysis）——檢定你抽取到的化合物是否比隨機更常落在某個路徑／模組中。 |
| InChIKey | 一種雜湊化的化學結構識別碼，用來辨識兩個 MS-DIAL 列是否為同一化合物。 |
| KEGG pathway | 各物種的 KEGG 反應目錄，將相關化合物分組。 |
| KEGG module | 各 Group 的 KEGG 目錄，描述化合物中更緊密的功能單元。 |
| adduct | 中性分子的某種特定離子化形式（例如 `[M+H]+`、`[M+Na]+`）。 |
| isotope peak | 化合物單一同位素（M0）峰的天然豐度 M+1 / M+2 衛星峰。 |
| hypergeometric test | ORA 背後的統計量——不放回抽取時恰好抽到這麼多路徑化合物的機率。 |
| universe / foreground | 完整候選化合物集合（母體，`N`）相對於受檢定的被抽取子集（前景，`K`）。 |
| volcano plot | 第二階段以 `log2(FC)` 對 `−log10(p_adjusted)` 繪製的散布圖。 |
| dot plot | 第三階段富集分析的結果圖。 |
| fold enrichment | 在 ORA 中某路徑相對於隨機被代表的倍數。 |
| PQN | 機率商數正規化（Probabilistic Quotient Normalization）——一種 NMR 風格、針對樣本稀釋的校正。 |
| arcsinh | 反雙曲正弦，一種變異數穩定化轉換，對大值表現得像對數，在接近零處則近似線性。 |
| measurable metabolome | 此檢測原則上能偵測到的化合物集合——富集分析母體的基礎。 |
| biosample | 一個生物樣本標籤，在雙模式中將正離子模式與負離子模式的進樣配對起來。 |

> **注意：核心概念 — 缺失值與真正的零。** 空白或 `"NA"` 儲存格會變成內部的「缺失」標記（`f64::NAN`），下游統計會**跳過**它；而寫成 `0` 的則是真正的零，會**參與**運算。詳見完整的 [缺失值 vs 真正的零](#missing-values-nan-vs-true-zeros-00) 章節。

  - [如何閱讀本手冊](#how-to-read-this-manual)
  - [管線一覽](#pipeline-at-a-glance)
  - [關鍵術語與符號](#key-terms--symbols)
  - [第一階段 — 輸入解析](#stage-1--input-parsing)
    - [MS-DIAL `.txt`](#ms-dial-txt)
    - [群組對應 `.csv`](#group-mapping-csv)
    - [第一階段 → 第二階段的關卡](#stage-1--stage-2-gate)
  - [第二階段 — 正規化、去除重複與 DAM](#stage-2--normalization-deduplication--dam)
    - [差異累積代謝物（DAM）](#differentially-accumulated-metabolites-dam)
      - [1. 未知特徵過濾（預設開啟）](#1-unknown-feature-filter-default-on)
      - [2. 逐特徵前置過濾](#2-per-feature-pre-filter)
      - [3a. 方法：Student t 檢定（等變異）\[參數法，**預設**\]](#3a-method-students-t-test-equal-variances-parametric-default)
      - [3b. 方法：Welch t 檢定（不等變異）\[替代參數法\]](#3b-method-welchs-t-test-unequal-variances-alternative-parametric)
      - [3c. 方法：Brunner–Munzel 檢定 + Cliff's δ \[無母數\]](#3c-method-brunnermunzel-test--cliffs-δ-non-parametric)
      - [4. 多重檢定校正](#4-multiple-testing-correction)
      - [5. 趨勢分類](#5-trend-classification)
      - [6. 火山圖](#6-volcano-plot)
      - [7. 將圖匯出為 PNG](#7-exporting-the-figure-as-png)
      - [DAM 值得注意的事項](#dam-caveats-worth-knowing)
    - [以 InChIKey 去除重複](#deduplication-by-inchikey)
      - [級聯決策表](#cascade-decision-table)
      - [稽核 CSV](#audit-csv)
      - [停用](#opt-out)
    - [樣本正規化](#sample-normalization)
      - [為何 Sum / Median / Metadata 重新縮放到中位數因子（而非除到固定常數）](#why-sum--median--metadata-rescale-to-the-median-factor-rather-than-divide-to-a-constant)
      - [生命週期](#lifecycle)
      - [啟動時的錯誤](#errors-at-startup)
      - [值得注意的事項](#caveats-worth-knowing)
  - [第三階段 — 富集分析（過度代表分析）](#stage-3--enrichment-over-representation-analysis)
    - [富集分析設定畫面](#enrichment-analysis-setup-screen)
    - [富集分析結果畫面](#enrichment-analysis-result-screen)
    - [路徑模式](#pathway-mode)
    - [模組模式](#module-mode)
    - [開始新一輪分析](#starting-a-new-analysis-round)
  - [進階主題與參考](#advanced-topics--reference)
  - [缺失值（`NaN`）與真正的零（`0.0`）](#missing-values-nan-vs-true-zeros-00)
  - [雙模式（正 + 負離子）輸入](#dual-mode-positive--negative-ionization-input)
    - [何時使用雙模式](#when-to-use-dual-mode)
    - [準備輸入](#preparing-inputs)
    - [不平衡或缺少某模式的樣本](#unbalanced-or-missing-mode-samples)
    - [第一階段 UI](#stage-1-ui)
    - [第二階段（共用設定、逐模式 DAM）](#stage-2-shared-setup-per-mode-dam)
    - [第三階段 — 雙模式 N 與 K 的運算](#stage-3--dual-mode-n-and-k-math)
    - [計算範例](#worked-example)
  - [快取與來源](#caches-and-provenance)
  - [儲存與載入工作階段設定（再現性）](#saving-and-loading-session-settings-reproducibility)
    - [檔案內容](#whats-in-the-file)
    - [完整檔案，逐欄位](#the-full-file-field-by-field)
    - [各按鈕何時可用](#when-is-each-button-available)
    - [載入流程](#loading-workflow)
    - [若我在上傳中繼資料前就載入設定，會怎樣？](#what-if-i-load-settings-before-uploading-metadata)
    - [手動編輯 JSON](#hand-editing-the-json)
  - [回報問題](#reporting-bugs)
  - [主要參考文獻](#key-references)

---

<a id="stage-1--input-parsing"></a>
## 第一階段 — 輸入解析（Stage 1 — Input parsing）

**平實說明：** 第一階段是你載入資料的地方。你把 MS-DIAL 匯出檔與一份小巧的群組對應試算表拖進來，metabolopan 會把它們解析成一張工作表格——並小心地讓*缺失*儲存格與*真正的零*保持區分。當兩個檔案都乾淨地解析、且你的群組有效時，**Continue to DAM** 按鈕就會亮起，你便進入第二階段。

metabolopan 以兩種模式之一接收輸入，依你載入的 MS-DIAL `.txt` 檔數量決定。
**單模式**是一個 MS-DIAL `.txt` + 一個群組對應 `.csv`；**雙模式**是兩個 MS-DIAL `.txt` 檔（一個正離子、一個負離子）+ 一個含 `biosample` 欄的群組對應 `.csv`。
下方的檔案格式對兩種模式皆適用；雙模式特有的機制（biosample 配對、群組對等性檢查、逐模式 DAM）涵蓋於下方的 [雙模式（正 + 負離子）輸入](#dual-mode-positive--negative-ionization-input)。

<a id="ms-dial-txt"></a>
### MS-DIAL `.txt`

- 前 4 列是 MS-DIAL 的中繼資料（`Class`、`File type`、`Injection order`、`Batch ID`）；第 5 列是欄位標題。
  當某一欄的 `File type` 值非空白、且不是 `"NA"`、也不是字面上的列標籤 `"File type"` 時，該欄會被視為真正的樣本進樣——並保留在 `sample_cols` 中。
  這**包含** `Sample` 與 `Blank`（製程空白）；只排除 MS-DIAL 各群組的 `Average` / `Stdev` 彙總欄（標記為 `NA`）。
- **版本相容性。** 同時支援 MS-DIAL 4 與 MS-DIAL 5 的 Alignment 匯出檔。
  欄位是依名稱查找，因此 MS-DIAL 5 重新排序／改名的評分欄（它把 `Dot product` 拆成 `Simple` / `Weighted dot product`）也能以相同方式解析；metabolopan 只使用兩個版本共有的欄位。
- **缺失值。** 空白／僅含空白字元／`"null"`／`"NA"`／無法解析的強度儲存格會變成 `f64::NAN`——這是內部的「缺失」標記。
  明確寫成 `"0"` 的則維持 `0.0`。
  這可避免下游統計把缺測值與真正的零混為一談。
  詳見 [缺失值 vs 真正的零](#missing-values-nan-vs-true-zeros-00)。

<a id="group-mapping-csv"></a>
### 群組對應 `.csv`（Group mapping）

- CSV 必須包含一個名為 `sample` 的欄位與一個名為 `group` 的欄位，**位置與順序不限**。
  選用的 `biosample` 欄（位置不限；雙模式所必需——詳見下方 [雙模式輸入](#dual-mode-positive--negative-ionization-input)）會以名稱辨識。
  其後的任何欄位都會被解析為選用的中繼資料。
  欄名採精確比對（大小寫敏感）；缺少 `sample`/`group` 欄、重複的 `sample`/`group`/`biosample` 欄、空白的 `group` 儲存格、或重複的 `sample` 名稱，都會以明確的錯誤訊息拒絕。
  中繼資料欄會在載入時逐欄分類：若某欄的非空白儲存格全部都能解析為數字，就會出現在第二階段的 **Metadata column** 正規化單選按鈕中（例如 `dry_weight`、`dilution`、`total_protein`）；若某欄有任何非空白且非數字的儲存格（例如像 `CTR-01` 這樣的 `biosample` 標籤），則會被靜默地排除在該單選按鈕之外，並在應用程式內的日誌窗格以一行 WARN 說明被略過的是哪一欄、以及有多少儲存格無法解析。
  空白的中繼資料儲存格會解析為 `None`。
  出現在 MS-DIAL `.txt` 中、卻不在 CSV 裡的樣本會被標記為 `Unassigned`；CSV 中指名了 `.txt` 所缺樣本的列則會記錄為警告並忽略，並在輸入頁以紅色橫幅列出。
  **未指派樣本只在第一階段可見**——輸入摘要面板會以黃色的 `Unassigned (N samples)` 列顯示它們，讓你知道它們存在，但當你在第二階段設定畫面按下 **Start DAM** 時，它們就會從工作矩陣中被捨棄。
  正規化、去除重複、DAM 統計或任何下游匯出都不會看到它們。
  若要讓某個樣本納入分析，請在中繼資料 CSV 中為它加上真正的群組標籤；若要完全排除某個樣本（連第一階段都不顯示），請從 MS-DIAL `.txt` 的 File type 列移除它那一欄（把該欄設為 `NA`）。

<a id="stage-1--stage-2-gate"></a>
### 第一階段 → 第二階段的關卡（Stage 1 → Stage 2 gate）

**Continue to DAM** 按鈕會維持停用，直到：

- 兩個檔案都成功解析；
- 第 #1 槽位的離子化模式單選按鈕已設定；
- 存在 ≥ 2 個不同的非 `Unassigned` 群組；
- 每個可指派群組都有 ≥ 2 個樣本（下游統計所必需）。

---

<a id="stage-2--normalization-deduplication--dam"></a>
## 第二階段 — 正規化、去除重複與 DAM（Stage 2 — Normalization, Deduplication & DAM）

**平實說明：** 以下全部都是單一個第二階段設定畫面上的選項；它們合起來設定一次 DAM 執行。
下方按照最重要的順序呈現——先是決定哪些代謝物有差異的核心統計檢定（DAM），接著是清理重複化合物列的去除重複，最後是在這一切之前校正技術性負載的選用樣本正規化。

<a id="differentially-accumulated-metabolites-dam"></a>
### 差異累積代謝物（DAM, Differentially Accumulated Metabolites）

**平實說明：** DAM 是第二階段的核心——它會在你所選的兩群組（分子組 vs 分母組）之間，逐一檢定每一個代謝物，告訴你哪些真的有差異。你挑一種統計方法（**Student**、**Welch** 或 **Brunner–Munzel**）、挑一種 FDR 校正（**BH** 或 **BY**）、點擊 **Start DAM**，就會得到一張火山圖，外加一張可匯出的「上調／下調／不顯著」特徵表格。

每個特徵都會在使用者所選的分子組與分母組之間被獨立檢定。
本軟體提供三種統計方法；它們都遵循相同的整體流程。

**我該挑哪種檢定？**（「離散程度」／變異數 = 一個群組的值有多分散。）

| Method | 在以下情況選它… | 性質 |
|--------|---------------|--------|
| **Student's t-test** *(預設)* | 兩群組離散程度相近、樣本數也相近 | 參數法 |
| **Welch's t-test** | 其中一組明顯比另一組更分散 | 參數法 |
| **Brunner–Munzel + Cliff's δ** | 即使經過 `arcsinh`，資料仍偏斜或呈現存在／缺失 | 無母數 |

<a id="1-unknown-feature-filter-default-on"></a>
#### 1. 未知特徵過濾（預設開啟）（Unknown-feature filter）

`InChIKey` 為 `null` 的特徵（MS-DIAL 的「Unknown」鑑定）會在任何統計工作之前被捨棄，這樣 FDR 校正的 `m` 就不會納入那些終究無法進入第三階段 ORA 的代謝物。
若使用者特別想對未鑑定特徵取得統計結果（例如標記出供後續鑑定的候選），可在第二階段設定中取消勾選 **Drop unknown features (no InChIKey)** 核取方塊。

<a id="2-per-feature-pre-filter"></a>
#### 2. 逐特徵前置過濾（Per-feature pre-filter）

對每個剩下的特徵，會先捨棄合併 `numerator ∪ denominator` 欄中的 NaN 值，然後依序要求：(i) 分子組有 ≥ 2 個非 NaN 值、(ii) 分母組有 ≥ 2 個非 NaN 值、(iii) 合併後的 `nunique > 1`、以及 (iv) 合併後的 `IQR > 0`。
未通過任一檢查的特徵會從結果中移除，並計入 UI 中可見的 `skipped` 計數。

<a id="3a-method-students-t-test-equal-variances-parametric-default"></a>
#### 3a. 方法：Student t 檢定（等變異）[參數法，**預設**]（Student's t-test）

古典（同質變異數）形式。
當各組樣本數相近、且兩組離散程度大致相當時最適用——在這些假設下，它比 Welch 略具檢定力。
**新工作階段的預設值**：搭配 **Log transformation**（arcsinh）步驟（同樣預設開啟），它是本專案的標準起點。
若你懷疑變異數不等，請改用 Welch；若分布偏斜到連轉換都不足以應付，請改用 Brunner–Munzel。

- 與 Welch 共用的選用前置檢定轉換：當第二階段設定的 **Log transformation** 核取方塊被勾選時（預設開啟；`SessionSettings.log_transform = true`），會對每個非 NaN 儲存格套用 `arcsinh(x)` 作為變異數穩定化步驟（asinh 能處理零／負值，而 log10 會把它們變成 NaN）。
  未勾選時，此步驟會被略過，原始工作矩陣的值會直接流入 t 檢定。
- 合併變異數 `sp² = ((na − 1)·va + (nb − 1)·vb) / (na + nb − 2)`，固定自由度 `df = na + nb − 2`，雙尾 p 值經由 Student-*t* CDF 求得。
- **倍數變化（FC）的尺度取決於 `log_transform`。** 原因：在 `log_transform=true` 之下，*t* 統計量是在 arcsinh 轉換後的尺度上計算，但 arcsinh 對正值是凹函數，因此由 Jensen 不等式可知，兩個重尾群組的*原始*平均比值，可能在**正負號**上與 *t* 檢定實際評估的 arcsinh 平均差不一致。
  若把原始平均比值與 arcsinh 尺度的 *p* 值並列回報，會靜默地誤判由離群值驅動的特徵（例如 `num=[0.1]×9 + [100]` vs `den=[5]×10`，得原始 FC ≈ 2.02 ⇒「Up」，但 Welch *t* ≈ −3.25、*p* ≈ 0.01 ⇒「Down」）。
  參數方法分支會讓尺度一致：
    - `log_transform=false`（原始尺度）：`FC = mean(numerator) / mean(denominator)`，`log2(FC) = log2(FC)`。
    - `log_transform=true`（arcsinh 尺度）：`log2(FC) = (mean(arcsinh(num)) − mean(arcsinh(den))) / ln(2)`，且 `FC = 2^log2(FC)`。
      在相同資料上，`log2(FC)` 的正負號**保證**與 *t* 統計量的正負號一致。
      對大的 *x*，arcsinh(x) ≈ ln(2x)，因此 `log2(FC)` 會漸近於 `log2(GM(num) / GM(den))`——即 limma / DESeq2 的古典對數倍數變化。
      對小的 *x*（接近 0），arcsinh(x) ≈ x，因此 `log2(FC)` 會退化為縮放後的算術平均差，而非真正的比值。
      這是變異數穩定化已被記載的結果；等價的對數 FC 詮釋只在大 *x* 漸近區（arcsinh 與 ln 對齊處）成立。
      CSV 匯出會經由 `fc_basis` 欄（`mean` / `median` / `arcsinh-mean`）標示目前作用的基準，讓下游使用者無需重跑流程即可辨識某個數字位於哪種尺度上。

<a id="3b-method-welchs-t-test-unequal-variances-alternative-parametric"></a>
#### 3b. 方法：Welch t 檢定（不等變異）[替代參數法]（Welch's t-test）

與 Student 屬於同一參數族，但不假設變異數相等。
當各組離散程度明顯不同時、或當你不確定而想用較安全的預設時，請用它。

- 與 Student 相同的選用前置檢定轉換（僅 `arcsinh`，由第二階段 **Log transformation** 核取方塊控制；預設開啟）。
- Welch 的 t 統計量是用（可選地經 arcsinh 轉換的）值、以 NaN 感知的平均數與變異數計算，搭配 Welch–Satterthwaite 自由度，再經由 Student-*t* CDF 轉換為雙尾 p 值。
- **倍數變化的尺度與檢定尺度一致**——規則與上方的 Student 相同。
  在 `log_transform=true` 之下，`FC` 位於 arcsinh 尺度，因此其正負號永遠與 Welch *t* 的正負號一致。
  在 `log_transform=false` 之下，`FC` 是古典的原始平均比值。

> **⚠ 警告：Welch / Student 的邊角情況 — 某組變異數為零。** 當某一組的每個重複樣本都是相同的值時（例如某特徵在某條件的所有樣本中都低於偵測極限、而被填補成一個常數），Welch–Satterthwaite 自由度會塌縮為*另一*組的 `n − 1`。
> 對 `n = 2` 而言，這會得到 `df = 1`，使 *t* 分布極寬、p 值極為保守——即使兩群組肉眼可見地分得很開時也是如此。
> 這是標準的數學行為（與 R 的 `t.test(var.equal=FALSE)` 和 SciPy 的 `ttest_ind(equal_var=False)` 完全一致），但對代謝體學而言，受影響的特徵往往對應到你可能想保留的、真正的「在某條件存在、在另一條件缺失」訊號。

`run_dam` 每次執行會發出一行 INFO 日誌，回報觸發此路徑的特徵數量（在你的工作階段日誌 `<data_dir>/metabolopan/logs/session_*.log` 中尋找 `zero_variance_features=N`，其中 `<data_dir>` 為 `dirs::data_dir()`——macOS `~/Library/Application Support`、Linux `~/.local/share`、Windows `%APPDATA%`）；當 N > 0 時，可考慮改用 Brunner–Munzel 方法重跑，它是以秩為基礎，會以不同方式處理這個邊角情況。

> **給開發者：** 這個診斷計數器使用相對容差——變異數低於
> `(max(|mean|, 1))² × 1e−20` 者即被標記——因此某組在浮點數雜訊範圍內為常數的特徵（例如在高強度尺度下對位元相等的正規化前輸入做算術，使得 `var ≈ ε² × c²` 非零、但自由度的病態行為仍以相同方式發作）也會計入。
> t 檢定函式本身內部的逐方法 `var == 0.0` 守衛維持不變——只有這個診斷計數器被放寬了。

<a id="3c-method-brunnermunzel-test--cliffs-δ-non-parametric"></a>
#### 3c. 方法：Brunner–Munzel 檢定 + Cliff's δ [無母數]（Brunner–Munzel + Cliff's δ）

當各組的強度分布偏斜或不相等、且變異數穩定化轉換仍不足以應付時最適用。
代謝體學資料常常難以用高斯假設妥善描述（高度偏斜的對數分布、頻繁的存在／缺失模式、批次假影），因此在這些情況下，Brunner–Munzel + Cliff's δ 能在工作流程所見的各種離散情形中提供更誠實的 p 值。
當預設的 Student t 檢定（即使經過 `arcsinh`）擬合不佳時——例如高度偏斜或以存在／缺失為主的特徵，或當你要對應先前已發表的無母數分析時——請透過第二階段設定的單選按鈕選用它。

- Brunner–Munzel 統計量是用 `numerator ∪ denominator` 上的中位秩計算，搭配類 Welch–Satterthwaite 的自由度，再經由 Student-*t* 分布轉換為雙尾 p 值。
  行為與 SciPy 的 `brunnermunzel(distribution='t')` 和 R 的 `lawstat::brunner.munzel.test` 一致——`sqrt` 內的 W 分母為 `nx·Sx + ny·Sy`。
- Cliff's δ 效應量：`(gt − lt) / (n · m)`，其中 `gt` 與 `lt` 分別是嚴格大於與嚴格小於的配對計數。範圍 `−1` 到 `+1`；|δ| ≥ 0.33 是此處採用的慣用「中等效應」門檻。

  > **範例：** Cliff's δ 是某一組隨機取一個重複樣本，其量測值高於另一組隨機取一個重複樣本的頻率。`δ = 0` 表示兩群組完全重疊；`|δ| = 1` 表示完全分離；`|δ| ≥ 0.33` 表示兩者約有三分之二的時候不一致（Cliff 1993）。

- 倍數變化使用群組**中位數**：`FC = median(numerator) / median(denominator)`，`log2(FC) = log2(FC)`。
  中位數對離群值穩健，與以秩為基礎的檢定哲學一致。

<a id="4-multiple-testing-correction"></a>
#### 4. 多重檢定校正（FDR）（Multiple-testing correction）

**平實說明：** 當你一次檢定數千個代謝物時，有些會純粹靠運氣看起來「顯著」——在 p < 0.05 下檢定 5,000 個特徵，即使什麼都不是真的，你也會預期約有 ~250 個假陽性。FDR 校正會把每個原始 p 值膨脹成一個考慮整個檢定家族的 q 值，以此加以節制。可以把 **BH** vs **BY** 單選按鈕想成一個敏感度／安全性轉盤：BH 是標準、發現較多的設定；BY 是較嚴格、在特徵彼此相關時更安全的設定。

每次第二階段執行都會對逐特徵的 p 值套用使用者所選的偽發現率（FDR）校正，不論這些 p 值是由哪種統計方法產生。
第二階段設定畫面提供一個有兩個選項的單選按鈕：

- **Benjamini–Hochberg (BH) procedure**——預設值。
  BH 假設檢定之間相互獨立或具正向迴歸相依（亦即假設這些檢定不會串通），並產生較多的發現。
- **Benjamini–Yekutieli (BY) procedure**——需主動選擇，較嚴格。
  把 BH 的 q 值乘上精確的調和因子 $c(m) = \sum_{i=1}^{m} \frac{1}{i}$（對大的 m 而言 ≈ ln(m) + γ，因此在 m = 5,000 時 BY 大約比 BH 保守 9 倍——那 ~9 倍就是你要付出的代價）。
  BY 在任意正向相依下都能控制 FDR，因此當許多特徵在生物上彼此相關時（例如共享路徑成員的代謝物），它是較安全的選擇。

NaN 的 p 值在任一方法下都會原樣以 NaN 通過校正。
所選的方法會回報在火山圖的註解條上（例如 `FDR(BH)<0.05`），並寫成每個 DAM CSV 匯出檔開頭的 `# FDR: BH` / `# FDR: BY` 註解行，因此螢幕截圖與下載檔都能自我說明。
參考文獻：Benjamini & Hochberg (1995)；Benjamini & Yekutieli (2001)。

<a id="5-trend-classification"></a>
#### 5. 趨勢分類（Trend classification）

它會在使用者調整門檻時即時重新計算——絕不儲存在結果中。
預設門檻：`FC = 2.0`（等同 |log2(FC)| > 1.0）、`FDR = 0.05`、`|δ| ≥ 0.33`（僅 BM）。

- Student / Welch（皆為參數法，無效應量）：`Up` 當且僅當 `FDR < threshold` 且 `log2(FC) > log2(fc_threshold)`；`Down` 當且僅當 `FDR < threshold` 且 `log2(FC) < −log2(fc_threshold)`。
  δ 門檻在參數檢定中被忽略。
- Brunner–Munzel：參數規則**且** `|δ| ≥ delta_threshold`。`δ = None` 的特徵（BM 因某組少於 2 個非 NaN 值而無法計算效應量）會被分類為 `NotSignificant`。

<a id="6-volcano-plot"></a>
#### 6. 火山圖（Volcano plot）

X 軸 = `log2(FC)`，Y 軸 = `−log10(p_adjusted)`。
**X 軸所代表的內容取決於作用中的方法與 `Log transformation` 切換**——對 `log_transform=false` 的 Welch / Student 是平均比值，對 `log_transform=true` 的 Welch / Student 是 arcsinh 平均差（以 log2 為單位），對 Brunner–Munzel 是中位數比值。
詳見上方第 3 節；作用中的基準會記錄在每個 `DamFeature` 上的 `fc_basis`（`mean` / `arcsinh-mean` / `median`）。
三種顏色依趨勢分類而定（紅／藍／灰，透明度 α ≈ 0.5）。
門檻線為實心黑色：位於 `−log10(FDR)` 的水平線、位於 `±log2(FC)` 的垂直線。
`log2(FC)` 為 `±∞` 的特徵（某組平均或中位數恰為 0）會被停靠在 X 軸邊緣 `±(xabs_max + 0.5)`，並加上小幅 jitter，使它們維持可見。
Y 軸上的對稱飽和：BH/BY q 值下溢到恰好為 `0.0` 的特徵（極大的 `|t|` / 極小的原始 p，在分得很開的群組中很常見）會被停靠在 Y 軸頂端（`y_max`）**正下方**，並加上每點向下、在 `−log10(q)` 單位下最多 `0.08` 的 jitter（與 X 軸 ±0.04 jitter 慣例尺度匹配），使多個飽和特徵不會堆在單一像素上。
底層的 `neg_log10_p_adjusted` 值仍為 `f64::INFINITY`，在 CSV 匯出中也照此記錄——只有螢幕上的位置被 jitter。
Y 軸在其他情況下僅為顯示用途而被裁切於 `finite_max + 1`；底層數值在 CSV 匯出中保留。
`NaN` 的 `neg_log10_p_adjusted` 保留給真正「p 無法計算」的情況（BM 下完全分層的群組；NaN-drop 後 `n < 2` 的參數檢定）——這些點會從圖中略去，但仍列在 CSV 裡。
X 軸標籤下方的單一註解條會摘要方法、作用中的 FC 基準（`FC: mean` / `FC: median` / `FC: arcsinh-mean`）、作用中的門檻，以及 ±∞ 計數——例如 `Method: Brunner-Munzel | FC: median | FDR(BH)<0.05, FC≥2, |δ|≥0.33 | −∞: 12  +∞: 8`。

**BM 點大小編碼 Cliff's δ 的量值。** 在 Brunner–Munzel 渲染中，每個散布點的半徑由該特徵的 `|Cliff's δ|` 對應而來：`|δ|=0` 給出仍可見的最小點，`|δ|=1` 給出約為預設半徑 1.3× 的點，中間量值則在兩個錨點之間線性縮放。
右側圖例會在既有的趨勢計數下方長出第二個 `|δ| size` 區塊，含三個位於 `|δ|=0/0.5/1.0` 的中性灰參考點——以這些參考點對散布點做大小比對，即可從圖上讀出量值。Welch / Student 渲染在全圖維持一致的點半徑，且**不**繪製 `|δ| size` 圖例區塊（那些檢定不產生可編碼的 Cliff's δ）。`|δ|` 未定義的 BM 特徵（某組非 NaN 值 `n < 2`）會退回預設半徑，並仍以適當的趨勢顏色渲染。

<a id="7-exporting-the-figure-as-png"></a>
#### 7. 將圖匯出為 PNG（Exporting the figure as PNG）

相同的三個匯出控制項位於預覽上方，此處與第三階段點圖畫面皆然：**Width (in)**、**Height (in)** 與 **DPI**。
它們精確地定義所儲存的影像：`pixel width = round(Width × DPI)`、`pixel height = round(Height × DPI)`（各自夾擠於 `[64, 20000]` px）。
`Width` / `Height` 範圍為 `1.0–40.0` in；`DPI` 範圍為 `72–1200`。第二階段預設值為 `3.5 × 2.2 in @ 300 DPI`（→ `1050 × 660` px）。

- **Width / Height（英吋）** 設定圖在頁面上的物理尺寸。`DPI` 值也會寫進 PNG 的 `pHYs` 區塊（每公尺像素數），因此版面工具（Word、InDesign、LaTeX `\includegraphics`）會把影像精確放置成那麼多英吋，而不是從原始像素數推斷尺寸。
- **DPI** 設定解析度：提高它會在*相同*物理尺寸下讓點陣更銳利（更多像素）——`300` 是線稿的常見期刊下限，`600` 則用於印刷品質。縮小 `Width` / `Height` 會讓圖變小；提高 `DPI` 則維持版面尺寸但增加細節。
- **一切一起縮放。** 字型、軸刻度、門檻線與散布點都相對於畫布調整大小，因此改變三者中的任何一個都會讓整張圖均勻地重新縮放——文字絕不會相對於圖變得過小或過大。（點圖只以 `Width × DPI` 來決定字型大小，因為它的高度會依列數自動配合——見*路徑模式*下的第 10 項。）

**所見即所得。** 螢幕上的預覽與下載的 PNG 來自*相同*的渲染器、以*相同*的尺寸產生——沒有另一道「匯出品質」處理。預覽*就是*那個檔案：相同的版面、字型、顏色、點位置與像素尺寸（你的螢幕可能會把它縮放以符合視窗，但儲存的像素與預覽的相符）。
預覽是你上一次 **Draw volcano** / **Re-draw volcano** 所產生的影像。改變門檻會把它清空（按鈕還原為 **Draw volcano**）；改變匯出尺寸則不會——因此在調整 **Width** / **Height** / **DPI** 之後，請點擊 **Re-draw volcano**，使預覽與 **Download volcano PNG** 將寫出的內容一致。

<a id="dam-caveats-worth-knowing"></a>
#### DAM 值得注意的事項（DAM caveats worth knowing）

- BM 以中位數為基礎的 FC 意味著小 n 研究（例如每組 3 個樣本）比 Welch 以平均數為基礎的 FC 更可能產生 `±∞` 的 log2(FC)，因為三個樣本中只要有一個零，就會把群組中位數拉到零。
  註解條會顯示 ±∞ 計數，因此這絕不會悄無聲息地發生。
- 每組 n = 2 的參數 t 檢定（Student 或 Welch）只有約 1–2 個自由度，並不可靠；第一階段關卡的「每組 ≥ 2 個樣本」要求讓你維持在地板之上，但要有充分檢定力的參數檢定會想要每組 ≥ 5 個。
  當等變異假設成立時，等樣本數的 Student 是三者中最敏感的；當變異數肉眼可見地不同時，Welch 是穩健的退路。
- 趨勢分類取決於作用中的門檻。
  火山圖與 CSV 匯出器都以相同的即時門檻分類每個特徵，因此剛繪製的圖與 CSV 在構造上一致。
  沒有自動重繪：改變任何門檻會把火山圖清空（按鈕還原為 **Draw volcano**），直到你重繪它，因此螢幕上的圖絕不會留下顯示過時的分類——它要嘛與目前的門檻相符，要嘛什麼都不顯示。

<a id="deduplication-by-inchikey"></a>
### 以 InChIKey 去除重複（Deduplication by InChIKey）

**平實說明：** MS-DIAL 經常把同一個化合物回報為好幾列（不同的加成物、同位素峰，或分裂峰）。去除重複會把每一組同化合物的列塌縮成單一個最佳特徵，使它們不會讓你的檢定家族倍增。它以第二階段設定畫面上的核取方塊形式執行（預設開啟），事後你可以下載一份稽核 CSV，看清楚究竟哪些列被捨棄、又是為什麼。

MS-DIAL 經常為解析到同一個化合物的多個 Alignment ID 各自輸出一列。
這有三種生物／儀器層面的成因：

1. **加成物多重性。** 同一個中性分子在正離子模式下會以 `[M+H]+`、`[M+Na]+`、`[M+NH4]+`、… 形式離子化（或在負離子模式下以 `[M-H]-`、`[M+Cl]-`、`[M+FA-H]-`、…）。
   每個加成物都會產生自己的 Alignment ID，但共用同一個 InChIKey。
2. **同位素峰。** MS-DIAL 會為 M0 單一同位素峰、以及 M+1 / M+2 天然豐度同位素峰各自輸出獨立的列（由 `Isotope tracking weight number`、或 `Adduct type` 中的 `[M+1]` / `[M+2]` 後綴標示）。
3. **層析峰分裂。** 當峰偵測不夠理想時，單一個高斯沖提峰可能被切成兩個相鄰的 Alignment ID，它們共用每一項鑑定資訊，只在 `Fill %` / `S/N average` 上有差異。

把所有重複都餵進 DAM，會讓 FDR（偽發現率）的家族大小相對於真實化合物數膨脹 2–5 倍，侵蝕統計檢定力。
第三階段 ORA 不受這種特定的數量膨脹影響：它的前景 `K`（被抽取化合物數）與母體 `N` 都是以 InChIKey 收斂的*唯一* KEGG 化合物集合，因此加成物、同位素峰與分裂峰這些共用 InChIKey 的重複都會塌縮成單一化合物，不論是否去除重複，`K` 都維持不變。
去除重複對第三階段仍然重要，但風險方向其實*相反*：某個低品質重複特徵的 DAM 趨勢若與其同伴相左，會使該共用化合物彙總為模式內 `Conflict` 而被*移出* `K`（是縮減前景、而非膨脹），而非讓級聯保留那個唯一的高可信特徵。

**去除重複以「預設啟用、可關閉」的切換開關呈現於第二階段設定畫面（預設開啟）。** 此級聯*純粹*是去除重複作業，並非通用的品質過濾器——`inchikey = None` 的特徵會原封不動地通過，而單一條目（一個 InChIKey 只對應一個 Alignment ID）即使鑑定品質不佳也會被保留。

<a id="cascade-decision-table"></a>
#### 級聯決策表（Cascade decision table）

在每個相同 InChIKey 的群組內，存活的特徵由此級聯中第一個能區分兩者的層級決定：

| Level | Field                              | Rule                                                                                                  |
|-------|------------------------------------|-------------------------------------------------------------------------------------------------------|
| 1a    | `MS/MS matched`                    | `True` > `False` > 空白                                                                              |
| 1b    | `Total score`                      | 數值大者勝出（廠商計算的加權綜合分數，涵蓋所有光譜相似度指標，含 dot product） |
| 2     | 加成物類別                         | `Primary` > `NonPrimary` > `Dimer` > `Isotope`；在 `Primary` 之中，`[M+H]+` / `[M-H]-` > `[M+Na]+` / `[M+NH4]+` / `[M+K]+` / `[M+Cl]-` |
| 3a    | `Fill %`                           | 數值大者勝出（各樣本峰覆蓋率）                                                                |
| 3b    | `S/N average`                      | 數值大者勝出                                                                                           |
| 4     | `Alignment ID`                     | 字典序較小者勝出（決定性的最終判定）                                             |

加成物分類是決定性的且區分大小寫：`Isotope` 由 `Isotope tracking weight number > 0`、或加成物字串中的 `[M+<n>]` 後綴偵測得出；`Dimer` 由開頭的倍率（`[2M+H]+`、`[3M-H]-`、…）偵測得出；`Primary` 是封閉的允許清單 `{[M+H]+, [M+Na]+, [M+NH4]+, [M+K]+, [M-H]-, [M+Cl]-}`；其餘一切（包含缺少加成物儲存格的情況）皆為 `NonPrimary`。

<a id="audit-csv"></a>
#### 稽核 CSV（Audit CSV）

當 DAM 執行是在啟用去除重複的情況下產生時，底部面板的 **Data** 分頁會在第二階段結果畫面上顯示一個 **Download dedup audit (CSV)** 按鈕（在任何富集分析之後的畫面則不顯示）。
CSV 格式：

```
# Deduplication audit — generated by metabolopan
# Total dropped: <N>; total kept: <M>; null-InChIKey passthrough: <K>
dropped_alignment_id,inchikey,winner_alignment_id,decided_at,loser_value,winner_value
```

`decided_at` 欄會告訴你每一次捨棄是由哪個級聯層級所決定（`MsmsMatched` / `TotalScore` / `AdductClass` / `FillPercent` / `SnAverage` / `Tiebreak`）；`loser_value` 與 `winner_value` 則承載決定性欄位在兩側各自的內容（若該側為 `None` 則留空）。
在雙模式執行中，檔案會包含每個模式各一份報告，以 `# Mode: POS` / `# Mode: NEG` 標頭行分隔。

<a id="opt-out"></a>
#### 停用（Opt-out）

在第二階段設定畫面取消勾選 **Deduplicate features by InChIKey** 即可停用。
在未勾選的情況下，DAM 執行與導入本功能之前的行為逐位元相同——每一個輸入列都會抵達前置過濾、FDR 的 `m` 等於前置過濾後的數量、且 DAM 結果上的 `dedup_report` 為 `None`。

<a id="sample-normalization"></a>
### 樣本正規化（Sample normalization）

**平實說明：** 樣本正規化是一個選用的第一步，它在任何統計之前校正樣本之間的技術性負載差異（進樣體積、稀釋、乾重）。預設為 **None**，它是安全的、什麼都不改變。挑 **Sum** 或 **Median** 以校正進樣負載；挑 **Metadata column** 以正規化到像乾重這樣的量測量；挑 **Quantile** 以用於同一基質、重複數充足的研究；或挑 **PQN** 以做 NMR 風格的稀釋校正。

在進行任何逐特徵統計之前，使用者可選擇一種*樣本軸*（逐欄）正規化，以校正樣本之間的技術性變異（進樣體積、稀釋、乾重、總離子流）。
每次 DAM 執行開始時，矩陣都會從最初解析得到的原始強度值（`intensity_raw`）重新正規化一次；`intensity_raw` 永遠不會被更動，因此切換方法是無損的。
預設為 `None`，會逐位元保留先前的行為。

除預設值外，另提供五種方法：

- **Sum。** 每個樣本的因子 = 該樣本所有非 NaN 強度的總和。
  輸出
  $$ x^{\prime}_{[i, j]} = x_{[i, j]} \times \frac{\underset{j}{median}(f_j)}{sum_j} $$
  乘上各樣本總和的中位數可保留整體量級，讓 Welch / Student 路徑中選用的 `arcsinh` 步驟（由第二階段 **Log transformation** 核取方塊控制；預設開啟）維持在有用的範圍內。詳見下方。
- **Median。** 形式相同，改以各樣本的 NaN 感知中位數作為因子。
- **Metadata column。** 使用者從中繼資料 CSV 解析出的選用數值欄中挑一欄（例如 `dry_weight`、`dilution`）。
  每個儲存格會除以該欄對應該樣本的值，再以所有樣本的中位數值重新縮放。
  該指定中繼資料欄資料不完整時的行為：
  - *缺失值（空白儲存格）：* 該樣本會從分析中被**捨棄**——該樣本欄的每個儲存格都會標記為 NaN，讓 DAM 的 NaN 感知機制將其排除在逐特徵統計之外。
    第二階段設定畫面會在使用者按下 **Start DAM** 之前，以一行黃色警告列出將被捨棄的樣本。
  - *非正值（零或負值）：* 會明確報錯並指出有問題的樣本與欄位。
    零／負的中繼資料是資料輸入問題、而非「缺失值」，因此立即失敗（fail fast）才是正確做法。
  - *非數字儲存格：* 在 CSV 載入時即解析，並在抵達第二階段之前就報錯。
  - *群組前置檢查：* 在進行任何正規化工作之前，執行器會檢查：在捨棄沒有值的樣本之後，所選的分子組與分母組是否仍各自至少保有 2 個樣本。
    若否，錯誤橫幅會指出失敗的群組、欄位、剩餘數量、以及所需的最小值（`2`）。
- **Quantile Normalization。** 強制讓每個樣本的分布對齊到一個共同的參考（各秩位在所有樣本間的平均）。
  在已排序位置 `[k, k+t)` 出現並列的項目，會被指派為這 `t` 個秩位上參考值的**平均（MEAN）**——`mean(reference[k..k+t])`。

  > **注意：** 這是對 Smyth 在 Bioconductor 支援討論串 #1569（2003，<https://support.bioconductor.org/p/1569/>）一句話的**字面**解讀——他說並列項目應取「對應之合併分位數（pooled quantiles）的平均」。
  > 廣為部署的標準實作——包含 Smyth 自己的 `limma::normalizeQuantiles(ties=TRUE)`（預設值）以及 Bolstad 的 `preprocessCore::normalize.quantiles`——則改以平均秩查表搭配線性內插來處理並列，也就是回傳並列「中間秩」位置上的參考值。
  > 兩種解讀僅在 `t == 2` 的並列、或參考在局部呈線性時才一致，而在參考曲度較大時、對 `t ≥ 3` 的並列就會分歧（這在以低於偵測極限值填補的代謝體學樣本底部很常見，例如下方的計算範例）。
  > 因此 metabolopan 的輸出在這種情況下會與 preprocessCore 和 limma **兩者**都不同——這是刻意的設計，在此載明，好讓你能在知情的前提下與標準工具比較。

  > **範例：** 參考為 `[1.5, 7.5, 52.5, 502.5, 55000]`，在已排序位置 1–3 有三項並列，在此會得到 mean(7.5, 52.5, 502.5) = **187.5**；`limma::normalizeQuantiles(ties=TRUE)` / `preprocessCore::normalize.quantiles` 則回傳 `reference[2]` = **52.5**。
  > ![quantile-normaliztion-in-r](./figure/quantile-normalization.png)

  這個分歧與各樣本是否擁有相同的非 NaN 數量無關——它純粹取決於 `t ≥ 3` 的並列如何對應到曲度較大的參考；至於各樣本非 NaN 數量是否相等，是另一個面向，於下方說明。
  - **各樣本非 NaN 數量不等。** 當樣本擁有不同數量的非 NaN 儲存格時（例如缺測情形不一致），參考會建立在大小為 `K = max(n_j)` 的共同分數秩格點上，並將每個樣本已排序的值線性內插到該格點上。

    > **給開發者：** 這與 limma 的 `(r − 1)/(n − 1) ∈ [0, 1]` 機制一致。
    > 它可避免「較長的樣本主導高秩」這個錯誤——過去一個僅有 3 個非 NaN 的樣本，其最大值會被對應到參考的第 60 百分位（其 5 個位置中的 `reference[2]`），而非參考的第 100 百分位。
    > 當所有樣本的非 NaN 數量同為 `K` 時，每個分數秩都會落在整數格點索引上，內插路徑會退化為直接查表，輸出與本次變更之前我們所發布、僅支援等長的版本逐位元相同。
    > NaN 儲存格維持 NaN。
- **Probabilistic Quotient Normalization (PQN)。** 一種 NMR 風格、針對樣本稀釋的校正：它假設大多數特徵不應改變，因此某樣本相對於參考光譜的*典型*逐特徵比值，可估計該樣本的稀釋因子，再將其除掉。
  Dieterle 2006：先在內部做總和正規化；從所選群集（預設為 `All samples`，亦可選擇限定於某個指定群組）建立逐特徵的參考光譜；對每個樣本，計算其逐特徵商數相對於參考的中位數（略過參考為零、NaN、或樣本值為 NaN 的特徵）；再除以該因子並重新縮放。
  未指派樣本永遠不會抵達此階段（它們在第一階段 → 第二階段邊界就被捨棄了，因此無論是參考群集或逐樣本因子迴圈都不會看到它們）。
  若某個*已指派*的樣本仍產生退化的商數中位數（NaN 或 0），PQN 會中止並列出有問題的名稱——請改用其他正規化方法，或從 MS-DIAL `.txt` 的 File type 列移除該樣本。
  分派器的 INFO 日誌行會顯示一個 `reference_features_used=N` 欄位，讓你看到該群集實際錨定了多少個特徵作為 PQN 參考（亦即 `median(cohort) > 0` 者）相對於總特徵數——在不重跑流程的情況下，這對診斷 QC 稀疏度很有用。

<a id="why-sum--median--metadata-rescale-to-the-median-factor-rather-than-divide-to-a-constant"></a>
#### 為何 Sum / Median / Metadata 重新縮放到中位數因子（而非除到固定常數）（Why … rescale to the median factor）

這三種方法共用同一個驅動機制。
對每個樣本欄 *j*，它會計算一個純量因子 `f_j`——欄總和（Sum）、欄的 NaN 感知中位數（Median）、或樣本的正值中繼資料（Metadata）——然後把每個有限儲存格改寫為

  $$ x^{\prime}_{[i, j]} = x_{[i, j]} \times \frac{M}{f_j}, \; where \; M = \underset{j}{median}(f_j) $$

`× M` 這一項是刻意的設計。
總和正規化的教科書形式是單純的 `x / f_j`（或對 CPM 式計數乘上 `× 10^6`），這會強制把每個樣本拉到*每單位*尺度：對 Sum 而言，每欄總和都會變成 1（比例，約 1e-5 – 1e-3）；對 Median 而言，每欄中位數都會變成 1。
我們改為乘回 `M`，即**各樣本因子的中位數**，使每欄的總和（或中位數）落在 `M`——*典型*樣本的原始量級——而非 1。
樣本之間的技術性負載（進樣體積、稀釋、乾重）仍被均衡；只有絕對強度尺度被保留下來。

- *為何重要——下游的 `arcsinh`。* 預設的 **Log transformation** 是 `arcsinh`，它只有在 *x* 夠大時才表現得像對數（`arcsinh(x) ≈ ln(2x)`）；對接近 0 的 *x* 而言，它基本上是**線性的**（`arcsinh(x) ≈ x`）。
  把資料除成比例會把整個工作矩陣推進那個近線性區，使 arcsinh 的變異數穩定化效果崩潰，並讓 t 檢定退化成在線性尺度上比較一堆極小的數。
  把數值維持在強度尺度（約 1e4 – 1e7）能讓 arcsinh 停留在它有用的類對數區——即「永遠不要接近 0」的目標。
  對**所有**儲存格乘上相同的常數 `M`，在 Brunner–Munzel 的中位數比值中會抵銷，但在 `arcsinh` 之下**不會**（它是非線性的），因此這個重新縮放正是為了保護 Student / Welch + `arcsinh` 這條路徑——也就是目前的預設。
- *為何取因子的中位數（而非平均數）。* 中位數較穩健——單一個負載特別高的樣本無法把目標尺度往上拉——而且它讓*典型*樣本成為錨點：該樣本的 `f_j ≈ M`，於是 `f_j / M ≈ 1` 幾乎不會改變它，而偏離尺度的樣本則往它靠攏。
- *數字範例。* 三個樣本的欄總和為 `6, 15, 24`，得 `M = median(6, 15, 24) = 15`；做完 `x / sum_j × 15` 後，每欄總和都變成 **15**（樣本 A 放大 ×2.5、C ×0.625、B 不變）——而非變成 1。
  若以各樣本中位數 `2, 20, 200` 做 Median 正規化，則 `M = 20`，每欄的中位數都會變成 **20**。
  Metadata 也相同，只是 `f_j` 改為所選欄的值，得到「在中位數乾重下的強度」。

總結而言，metabolopan 的「先除再重新縮放到中位數」搭配的是 `arcsinh`，使正規化與廣義對數轉換在數值上保持相容。
所選的 `M` 會在分派器 INFO 日誌中以 `scaling_to_median_factor=…` 回報。

<a id="lifecycle"></a>
#### 生命週期（Lifecycle）

正規化的選擇——以及其他每一個設定參數——會在整個工作階段的生命週期內、跨越每一次導覽轉換而保留。
退回上一階段絕不會丟失你的選擇；你只是回到上一個畫面，先前所有的選擇都原封不動。
（若你在第一階段重新挑選檔案，而先前的分子／分母組在新的中繼資料中已不存在，第二階段會卡住關卡，直到你重新選擇有效的群組。）第三階段沒有獨立的正規化步驟——第三階段富集分析看到的，就是那份（已正規化的）工作矩陣。

<a id="errors-at-startup"></a>
#### 啟動時的錯誤（Errors at startup）

正規化會在 DAM 的 tokio 任務生成之前同步執行，因此任何失敗（例如 `Sample 'A2' is missing a value in metadata column 'dry_weight'`）都會立即顯示在紅色橫幅上。
只有當工作矩陣為有限值且形狀正確時，DAM 任務才會啟動。

<a id="caveats-worth-knowing"></a>
#### 值得注意的事項（Caveats worth knowing）

- *Quantile* 假設各樣本的分布*理應*相同。
  對於同一基質、重複數充足的研究（例如細胞萃取物）這是合理的，但對於跨組織或跨生物體的比較——生物本質上在分布層級就有差異——則不成立。
- *PQN* 對大多數 NMR 式的稀釋變異很穩健。
  所選的參考群集很重要：當研究有一個乾淨的基線群組時，以它作為 PQN 參考，往往比 `All samples` 產生更清晰的生物訊號。
  **PQN 對樣本品質很嚴格**：若某樣本的逐特徵商數中位數為 `NaN`（沒有可用於對照參考的特徵）或 `0`（其非參考零特徵中有半數以上恰為 0——通常是稀疏／類空白的樣本），會以錯誤訊息列出有問題的樣本名稱。
  請從中繼資料 CSV 移除這些樣本，或改用較寬容的方法（None / Sum / Median / Metadata / Quantile）。
- *Metadata* 的值必須為嚴格正值——除法與量級保留步驟都假設正值。
  零與負值會報錯，而非靜默地通過。
- *Sum/Median* 會完全保留樣本內的特徵比值；它們是同一種轉換的「縮放到量級」版本。
  兩者的差別在穩健性：Sum 對每個樣本中少數高強度離群值敏感；Median 則忽略它們。
<a id="stage-3--enrichment-over-representation-analysis"></a>
## 第三階段 — 富集分析（過度代表分析，over-representation analysis）

第三階段接收你在第二階段找到的差異累積化合物，並提出單一個生物學問題：*它們是否比偶然預期的程度更集中於某條已知路徑（或模組）？*
舉例來說，若你的「up」化合物有一半全都屬於糖解作用，那麼該路徑就是*過度代表*的——一個值得回報的訊號。
每個模式（路徑 / 模組）執行的是完全相同的統計機制；它們只在檢定所針對的 KEGG 化合物集合目錄上有所不同。

第三階段接收第二階段的 DAM 結果，並提問：*「在我這份差異累積化合物清單中，哪些 KEGG 條目被過度代表？」*——這裡的「條目」指的是**[一條 KEGG 路徑](https://www.kegg.jp/kegg/pathway.html)**（路徑模式）或**[一個 KEGG 模組](https://www.kegg.jp/kegg/module.html)**（模組模式）。
兩種模式在超幾何檢定、使用者所選的 FDR（BH 或 BY）、以及可量測代謝體母體上共用完全相同的機制；它們只在 ORA 所操作的化合物集合目錄上有所不同。

<a id="how-over-representation-analysis-works-here"></a>
### 過度代表分析在此如何運作（How over-representation analysis works here）

把你能量測並對應到 KEGG 的每一個化合物想像成**罐子裡的一顆球**。
其中有些球是有顏色的——它們屬於路徑 *P*。
接著你伸手抓出一把：你的*差異累積*化合物。
超幾何檢定問的是，你抓出的這把球中有色球的數量是否*多於盲目運氣*所能預測的程度。

四個白話量驅動整個檢定（在公式之前先定義好，方便你閱讀）：

| Symbol | Plain meaning |
| --- | --- |
| **`N`** | **背景母體**——*在此平台上你所有能量測*且也成功對應到 KEGG 化合物的東西。罐子裡的所有球。 |
| **`K`** | **前景**——*你的差異累積化合物*（你抓出的那一把）。 |
| **`m_p`** | *整個罐子*中有多少顆球屬於路徑 *P*（有色球）。 |
| **`k_p`** | *你抓出的那些球*（`K`）中有多少顆屬於路徑 *P*——你實際抓到的有色球。 |

> **範例：**
> `N = 300` 個可量測化合物，你抓出 `K = 30` 個差異累積化合物，而某條路徑裡含有 `m_p = 10` 個化合物。
> 若你的 `K` 個命中是隨機散布的，機率預測只有 `30 × (10 / 300) = 1` 個會落在該路徑中。
> 你實際命中了 `k_p = 5` 個。
> 五對上期望值一，這是強烈的過度代表——一個很小的 *p* 值。

第三階段自帶一個獨立於第二階段選擇的 FDR 校正單選按鈕——**預設為 Benjamini–Yekutieli（BY）**，對路徑／模組 ORA 而言是較安全的選擇，因為條目本質上會共用化合物（許多化合物出現在多條路徑中），這違反了 BH 的獨立性假設。對重視跨工具再現性的使用者，Benjamini–Hochberg 仍可選用。
第三階段點圖的色標標題與註解列都會標明目前作用的方法（例如 `-log10(FDR (BY))` / `FDR: BY`），而富集分析 CSV 開頭的 `# FDR:` 註解行也會記錄此選擇供下游解析。

<a id="enrichment-analysis-setup-screen"></a>
### 富集分析設定畫面（Enrichment Analysis setup screen）

設定畫面是你選擇*要檢定什麼*以及*要多嚴格*的地方。
你選定一個模式、一個 KEGG 範圍（路徑用物種、模組用生物群組 Group）、一個方向過濾器，以及條目大小／FDR 旋鈕，然後按下 **Run Enrichment**。

第三階段設定畫面是使用者挑選以下項目的地方：

- **Analysis Mode**（Pathway / Module），透過單選切換鈕。
  兩種模式的選擇以及它們已取得的 KEGG 快取在工作階段的整個生命週期內共存——在模式間切換是即時的，且絕不會重新取得你已抓取過的資料。
- **KEGG 範圍。** 路徑模式顯示一個可搜尋的物種選擇器，內含預先積極載入的 KEGG 生物清單；模組模式則顯示下方 *模組模式* 所述的 Level + Group 選擇器。
  選定一個物種（或 Group）會在此畫面內就地觸發對應的 KEGG 取得——一個帶有標題說明的小進度條會串流逐路徑（或逐模組 + ETA）的進度，無須離開設定畫面。詳見下文。
- **Include DAM features with direction**（`Both` / `Up only` / `Down only`）。
- **Minimum number of compounds detected in a pathway/module**（即「最小條目大小」過濾器；預設 `1`，範圍 `[1, 20]`）。
  在建立 FDR 家族**之前**，捨棄其母體限定化合物數低於此門檻的路徑／模組——標準解釋見[路徑模式 step 5](#pathway-mode)。
- **FDR correction**（ORA 預設為 BY 程序；BH 程序可供跨工具再現性使用——見上文）。
- **`Run Enrichment`** 按鈕（取得作業進行中時停用；停用狀態的滑鼠懸停提示會說明是哪個取得作業擋住了按鈕）。

<a id="enrichment-analysis-result-screen"></a>
### 富集分析結果畫面（Enrichment Analysis result screen）

這三個控制項位於*結果*畫面上，讓你能在看過資料後迭代調整圖形，無須走回設定畫面。

- **Enrichment FDR threshold**（預設 `0.05`）。
- **Minimum hit count**（FDR 後的顯示過濾器；預設 `1`）。
  控制點圖顯示上限的 Top N 輸入欄位也位於此畫面，讓你可以在看過資料後再迭代調整，無須回到設定畫面。
- **Top N pathways**（預設 `20`）。

<a id="pathway-mode"></a>
### 路徑模式（Pathway mode）

路徑模式將你的化合物對逐物種的 KEGG 路徑目錄進行檢定。
下方的管線會解析每個特徵的身分、建立可量測母體、對每條路徑執行一次超幾何檢定、校正多重檢定，並繪製點圖。

管線如下：

1. **身分解析（[PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)）。** 對每個通過第二階段前置過濾的特徵（不只是 DAM 顯著的那些），透過 POST 至 `compound/inchikey/property/InChIKey/CSV` 將其 `InChIKey` 解析為一個或多個 PubChem CID。
   每批最多 200 個 InChIKey。
2. **KEGG 化合物轉換（[KEGG REST](https://www.kegg.jp/kegg/rest/keggapi.html)）。** 對每個唯一的 CID，透過 `/conv/compound/pubchem:CID1+CID2+...` 解析為一個 KEGG 化合物（`cpd:Cxxxxx`）。
   每批最多 10 個 CID，並受 KEGG 客戶端節流（每次請求間隔 334 ms，約 3 req/s，符合 KEGG 文件所載的軟性上限）。
   對應到 `glycan:` 或 `dr:` 的行會被濾除——只保留 `cpd:` 目標。
   HTTP 403 被視為速率限制訊號，並以 5 s 退避最多重試 5 次。
3. **多重對應規則。** 一個特徵就是一個化學實體。
   若 PubChem 對某個 InChIKey 回傳多個 CID（立體／位置／鹽類歧義）且它們全都解析到同一個 KEGG cpd，則該特徵對 DAM 化合物集合 `K` 與母體 `N` 各貢獻該 cpd **一次**。
   若它們解析到確實不同的 cpd，則每個 cpd 各自計數。
   若某特徵的 InChIKey 沒有 PubChem CID、或其 CID 全部無法對應到 KEGG cpd，則該特徵會從 `K` 與 `N` 中被丟棄，並顯示於底部面板 **Data** 分頁的對應漏斗中（`<N> InChIKeys → <N> PubChem CIDs → <N> KEGG cpds`）。
4. **母體定義（N）。** 母體是所有「通過第二階段前置過濾**且**成功經由 PubChem 與 KEGG conv 對應」之已註解特徵的唯一 cpd ID 的聯集——即此平台上的*可量測代謝體*。
   我們刻意採用「僅可量測」母體，使 p 值更能反映你的資料原本能說明的範圍。
5. **FDR 前的條目大小過濾。** 在任何超幾何工作之前，每條路徑的 `m_p` 會與使用者可調的 `min_entry_size`（預設 `1`，範圍 `[1, 20]`）比較。
   `m_p < min_entry_size` 的條目會從本次執行中**完全捨棄**——它們不對 FDR 家族貢獻任何 p 值、不出現在 CSV 中、也不出現在點圖上。
   被捨棄的數量會顯示於底部面板 **Data** 分頁的一行保留率資訊 `Tested: <surviving> / <total> (≥ N compounds in universe)`（模組模式為 `Tested: <surviving> (≥ N compounds in universe)`）。
   預設 `1` 讓前置過濾保持寬鬆——只有 `m_p = 0` 的條目會被捨棄（這種條目本來就只可能得到 `p = 1.0`），因此每條至少含一個可量測化合物的路徑都會被檢定。
   把它提高到 `3` 可額外排除 `m_p ∈ {1, 2}` 的條目，這些條目在典型的 `K`／`N` 值下在數學上無法被檢定——例如 `m_p = 1` 的條目最多只能產生 `k_p = 1`，得到原始 `p ≈ K/N`，這很少低於 `α = 0.05`，更少低於 BH 臨界值 `0.05/m`。
   這個取捨是對稱的：較低的 `min_entry_size` 會檢定更多路徑，但會擴大多重檢定家族 `m`。

   > **注意：** `m_p`（此處）與超幾何的 `m` 參數都使用交集的**集合**基數：在某個 KEGG 條目的 COMPOUND 區塊中列出多次的化合物只計**一次**，而非按出現次數計。
   > 這個 `min_entry_size` 旋鈕與 *Minimum hit count* **正交**：前者是**縮減 `m` 的 FDR 前條目過濾器**；*Minimum hit count* 則是**不改變 p 值的 FDR 後顯示過濾器**。

6. **逐路徑超幾何檢定。** 對每條通過條目大小過濾而存活的路徑 `p`，以
   `m_p = |unique(pathway.compounds) ∩ universe|`（該路徑落在可量測母體內之唯一 cpd ID 的集合基數——單一 COMPOUND 區塊內的重複 cpd ID 不會膨脹 `m_p`）以及
   `k_p = |K ∩ pathway.compounds|`：
   - `p_value = 1 - HypergeometricCDF(k_p - 1; N, m_p, K)`（觀測到*至少* `k_p` 個命中的上尾機率）
   - 若 `k_p, m_p, K, N` 中任一為零，實作會短路為 `p_value = 1.0`（避免未定義的 CDF 參數）。

   > **給開發者：** 實作還以 `debug_assert!` 強制 `K ⊆ N`（任何讓 K 洩漏出 N 之外化合物的上游退化都會在 dev/test 中大聲被捕捉；release 建置會在每次執行發出一條 `ERROR` 日誌，彙總任何超幾何定義域錯誤，使「所有條目皆不顯著卻無任何診斷」這種失效模式無法悄悄出貨）。

   - **Fold enrichment（富集倍數，效應量）。** 在每個 p 值之外，ORA 還會記錄一個效應量指標——觀測命中數除以虛無假設下的期望命中數。*期望值* = 若你的 `K` 個命中是隨機散布的，僅憑此路徑的大小有多少會落在其中：`expected_p = K · (m_p / N)`，故 `fold enrichment = k_p / expected_p = (k_p · N) / (m_p · K)`。
     `> 1` 代表過度代表（命中數多於機率預期）、`= 1` 恰如預期、`< 1` 代表代表不足。
     它是點圖的 **X 軸**，也是匯出 CSV 的 `Expected` / `EnrichmentRatio` 兩欄；它**只是效應量、不含顯著性**——一個只含單一化合物的條目（`m_p = 1`）靠一次幸運命中就能得到很大的 fold enrichment，這正是為什麼選取依 FDR 而非 fold enrichment（step 9）、以及為什麼存在 `min_entry_size` 前置過濾（step 5）的原因。
     邊緣情況：當 `expected_p = 0`（該條目中沒有可量測化合物）時，比值為未定義——內部為 `NaN`，在 CSV 中寫成**空白**儲存格。
7. **使用者所選的 FDR 校正**，透過第三階段設定畫面的獨立單選按鈕（預設為 Benjamini–Yekutieli 程序；Benjamini–Hochberg 程序只需一鍵切換；`None` 作為第三個選項，僅供探索性執行使用，見下文）。
   此單選按鈕刻意與第二階段的選擇相互獨立：兩個階段有不同的相依性樣態，許多使用者會合理地想要第二階段 BH（火山圖的跨工具再現性）+ 第三階段 BY（對共享化合物條目採保守 ORA）。
   對路徑／模組 ORA 我們**預設為 BY**：路徑大量共用化合物（糖解作用 ↔ TCA 共用 G6P、丙酮酸等），因此 BH 底層的獨立性假設被違反。多數生物學工具預設為 BH；若你需要可跨工具比較的 q 值，請把單選按鈕切到 BH。
   BY 在相依情況下較為保守；預期校正後的 p 值會一律較高（較不顯著）。
   `None` 完全跳過多重檢定校正——結果表與 CSV 中的 `fdr` 欄位會原封不動地攜帶原始 p 值。

   > **⚠ 警告：** 請**僅**將 `None` **用於探索性排序**，絕不可用於已發表的顯著性主張；在典型的 KEGG 路徑目錄上（約 300 條路徑受檢定），你純粹靠機率就會在 `p < 0.05` 預期約 15 個偽陽性。

   第二階段 DAM 單選按鈕**不**提供 `None`——約 13 k 個特徵的原始 p 會淹沒結果集；攜帶 `dam_fdr_method=NoCorrection` 的手工快照會防禦性地被強制改回 BH，並附帶一個 `tracing::warn!` 事件。

   **色階。** 每個標記的填色將 `-log10(FDR)`（`None` 時為原始 `-log10(p)`）編碼於 ColorBrewer **YlOrRd** 9 階漸層上——所顯示之最不顯著條目（位於顯示門檻的 FDR）為最淡的黃色，加深至最顯著者的深紅色；點與色標圖例共用單一個 `-log10` 跨距，因此相同顏色在兩者間代表相同的顯著性。
   目前作用的方法會記錄於點圖的色標標題（`-log10(FDR (BH))` / `-log10(FDR (BY))` / `None` 時為 `-log10(p-value)`——外層包裝被去除，因為軸上的值「就是」原始 p、而非 q），並記錄於匯出之富集分析 CSV 開頭的 `# FDR: BH` / `# FDR: BY` / `# FDR: None` 行。
   CSV 還攜帶額外的自我說明註解行，記錄該次執行所用的門檻：`# MinEntrySize: N`（FDR 前的條目大小過濾），以及在模組模式下 `# MinGroupOverlap: N`（Group 重疊門檻）。
   點圖本身在 X 軸下方還附帶一個四行的純文字註解區塊，讓審閱者僅憑圖即可重建 FDR 家族：

   ```
   Background universe = <N> compounds measured and mapped to KEGG
   Compounds of interest = <K> differentially abundant (increased | decreased | both directions)
   Pathways tested = <m>[ of <total>  ·  <dropped> skipped (< <min_entry> compounds each)][; ≥ <min_hit> hits required]
   Significance: FDR-adjusted, Benjamini–Yekutieli (BY)
   ```

   （當那些方法作用時，最後一行會讀作 `… Benjamini–Hochberg (BH)`，或 `raw p-value (no FDR correction)`。）
   `N` / `K` / `m` 等符號刻意完整拼出而非縮寫；受檢定數 `<m>` 是抵達 BH/BY 的條目數，也是每個原始 p 值被乘上的除數。
   `m` 分母等於**通過 FDR 前 `min_entry_size` 過濾**（step 5）的路徑數——即 `m = entries.len() − entries_dropped_by_min_entry_size`。
   協調器層級的 Group 過濾（模組模式）在更早的一層套用；到 FDR 執行時，`m` 已反映了這兩道過濾。
8. **顯示過濾（FDR 後）。** 一個使用者控制的 `min_hit_count`（預設 1）會把命中數較少的路徑從點圖與 CSV 中隱藏。
   這是一個*顯示*過濾器——`m` 已在所有存活條目上計算完畢，因此不論此設定為何，FDR 值都是誠實的。
   與 step 5 的 `min_entry_size` 不同：那一個是**縮減 `m` 的 FDR 前條目過濾器**；這一個是**不改變 p 值的 FDR 後顯示過濾器**。
9. **點圖的「選取」與「排序」——兩種不同的依據。** 點圖以**刻意不同的準則**選擇*要繪製哪些*條目以及*如何*把它們堆疊在 Y 軸上：
   - **選取（哪些條目出現）依統計顯著性。** 在通過 `fdr < threshold` 與 `min_hit_count` 過濾（step 7–8）的條目中，圖保留 **FDR 最低的 Top N**（`top_n`，預設 20，可在結果畫面調整）。
     因此所顯示的條目永遠是*最顯著*的那些——它們**絕不**依富集倍數選取。
   - **垂直順序（Y 軸）依效應量。** 被保留的條目接著依**富集倍數（觀測／期望）由大到小**排列，使富集倍數最大的條目位於**最上方**，整張圖沿 X 軸（X 軸本身即為富集倍數）讀起來像一道「大者在上」往下的階梯。
     平手時以 FDR 打破（較顯著者在前），再以條目 ID 打破。
     這符合 clusterProfiler 以 X 軸指標排序 Y 軸的慣例。

   > **注意：** 一條*極小*的路徑可能靠一次幸運命中就顯示出巨大的富集倍數——因此顯著性／FDR 決定**什麼**會出現，而富集倍數只是把存活者堆疊起來。
   > 請據此判讀點圖：**顏色與垂直位置 = 你有多確定**（顯著性）；**X 軸 = 效應有多強**（富集倍數）。

   實務上的後果：當顯著條目多於 `top_n` 時，被省略的是**最不顯著**（FDR 最高）的那些——*而非*富集倍數最小的。
   顯著性把關「是否納入」；效應量只負責排列已納入者。
   匯出的 CSV 與此無關：它列出每個存活條目並依 FDR 由小到大排序，附完整（未截斷）的名稱。
10. **點圖畫布高度。** 匯出圖的高度會自動配合實際顯示的列數——`clamp(min(top_n, displayed) × 0.3 + 1.0, 2.0, 40)` 英吋——並在你每次 Draw / Re-draw 時**重新計算**。
    因此若某次執行在你最初的 FDR 門檻下不顯著，而你在結果畫面放寬門檻並重繪，畫布會增大以容納新顯現的列，而非把它們塞進一張矮圖中（那會截斷 Y 軸標籤）。
    編輯 **Height (in)** 欄位會把它變成手動覆寫，並一直維持到下一次富集分析執行／重跑重置自動配合為止。

    **文字大小與條目數無關。** 標籤、軸標題、色標與 Hits 圖例會隨繪圖**寬度**（固定的 `Width (in) × DPI`）縮放，*而非*自動配合的高度——因此兩個條目的結果，其文字渲染大小與二十個條目的完全相同。
    高度 `2.0` 英吋的下限存在的原因，是讓全尺寸圖例在那些稀疏結果上總能容於畫布。
11. **匯出點圖（PNG 大小 + DPI）。** `Width (in)` / `Height (in)` / `DPI` 控制項與所見即所得保證的運作方式與火山圖完全相同——共用的 `pixels = round(inches × DPI)`、`pHYs` 物理尺寸、夾限、以及與預覽相同的渲染機制，描述於第二階段的[7. 將圖匯出為 PNG](#7-exporting-the-figure-as-png)。
    點圖特有的事實是：
    - 匯出預設為 `3.5 × 7.0 in @ 300 DPI`（`7.0` 是預設 `top_n = 20` 的自動配合高度）。
    - **Height** 自動配合所顯示的列數，並在每次 Draw / Re-draw 時重新計算，除非你覆寫它（上述第 10 項），而 `Width` 與 `DPI` 是你設定的普通數值。
    - 字型以 `Width × DPI` 為準，因此改變 `Width` 或 `DPI` 會重新縮放文字；改變 `Height` 則不會。

    預覽是你上次 `Draw dot plot` / `Re-draw dot plot` 所得的影像；在改變任何尺寸（或 `Top N`、FDR 門檻、最小命中過濾）之後，請點擊 `Re-draw dot plot`，使它與下載將產生的結果相符。

<a id="module-mode"></a>
### 模組模式（Module mode）

模組模式針對 KEGG *模組*（小型的功能性反應單元）而非整條路徑進行檢定，並以生物**群組 Group** 而非單一物種來界定範圍。
目錄選擇之後的所有下游——PubChem 對應、超幾何檢定、FDR、點圖——都與路徑模式完全相同。

模組模式執行與路徑模式完全相同的 PubChem → KEGG conv → 超幾何 → 使用者所選 FDR 流程，但 **(a)** 條目目錄是 KEGG 模組的集合、而非逐物種的路徑，且 **(b)** 使用者挑選的是一個**[生物群組（organism Group）](https://www.kegg.jp/kegg/tables/br08606.html)** 而非單一物種。
當某模組的 KEGG `COMPLETE` 區塊包含至少 `min_group_overlap`（預設 `1`）個來自所選 Group 的生物時，該模組就會被納入分析；這就是逐物種框架如何對應到全域模組目錄的方式。

1. **生物群組選擇。** 當 Analysis Mode 切換鈕設為 Module 時，第三階段 **Enrichment Analysis setup** 畫面會浮現一個 Level 單選按鈕（1 / 2 / 3）與一個 Group 下拉選單。
   在 Group 下拉選單正下方，一個 **Minimum group overlap** 數值控制項設定 `min_group_overlap` 門檻（預設 `1`，範圍 `1`–`min(Group size, 20)`）；其效果見下文的「模組 → Group 過濾」。
   Level 索引進 [KEGG 譜系欄位](https://www.kegg.jp/kegg/tables/br08606.html)（Level 1 為 `Eukaryotes`，Level 2 為 `Animals` / `Bacteria` 等，Level 3 為 `Mammals` / `Insects` 等）。
   KEGG 目前有 11,744 個生物，全都恰好有 4 個譜系層級；我們公開前三層。
   挑選一個 Group 會具現化 `org_codes`：屬於該 Group 的 KEGG 生物代碼集合（`hsa`、`ath`、…）。

2. **模組 → Group 過濾（[KEGG REST](https://www.kegg.jp/kegg/rest/keggapi.html)）。** 每個模組的 `/get/<module-id>` 回應攜帶一個 `COMPLETE` 區塊，列出該模組完整組裝的生物。
   當以下條件成立時，模組會被保留供 ORA 使用：
   ```
   |module.complete_orgs ∩ group_orgs|  ≥  min_group_overlap
   ```
   預設 `min_group_overlap = 1` 是寬鬆的（∃-重疊：Group 中任何單一生物就足夠）。
   較高的值會收緊過濾——例如 `min_group_overlap = 5` 要求 Group 的生物中至少有 5 個完整組裝了該模組。
   目前作用的門檻透過第三階段設定畫面上的 **Minimum group overlap** 控制項設定，並記錄於匯出 CSV 的 `# MinGroupOverlap:` 註解行，因此你發表的任何數字僅憑標頭 + 快取快照即可再現。

3. **母體與 K——同路徑模式。** PubChem 與 KEGG-conv 階段與模式無關。
   `N` 仍是可量測代謝體（成功對應到 KEGG cpd 的 DAM 特徵）；`K` 仍是符合目前作用之方向過濾（Up / Down / Both）之 DAM 特徵的 cpd 集合。
   模組模式*不會*以「所有模組化合物」或「所有 KEGG 化合物」替代 `N`。

4. **逐模組超幾何檢定。** 與路徑模式相同：對每個被保留的模組 `m`，
   `M_m = |module.compounds ∩ universe|`、`k_m = |K ∩ module.compounds|`，以及
   `p_value = 1 - HypergeometricCDF(k_m - 1; N, M_m, K)`，並使用相同的零輸入短路。

5. **使用者所選的 FDR 校正**——選項與預設皆同路徑模式（對共享化合物條目預設 BY；BH 可選）。
   `m` 分母等於**被保留模組**（Group 過濾之後）的數量，而非 KEGG 目錄中約 573 個模組的總數。
   這是正確的虛無假設：ORA 問的是「在*可能*適用於此生物群組的模組中，哪些被過度代表？」把分類上無關的模組納入 `m` 會在不貢獻生物學訊號的情況下把 FDR 往上扭曲。

6. **空 COMPOUND 模組計數器。** 有些 KEGG 模組（signature／reaction-only 模組）根本沒有 `COMPOUND` 區塊。
   當 `compounds = []` 時該模組的 `M_m = 0`，因此——就和路徑模式中任何 `M_p = 0` 的條目一樣——它會在任何超幾何檢定之前被 FDR 前的 `min_entry_size` 過濾器捨棄：它永遠不會走到 `p_value = 1.0` 的短路，也不對 FDR 家族貢獻任何 p 值。
   一個獨立的空 COMPOUND 計數器仍會統計它們，底部面板 **Data** 分頁會以一行 `With compound list: <kept>  (−<empty> empty)` 呈現，使默默捨棄絕不侵蝕信任。
   （對等的路徑模式回報已列入規劃。）

**模組模式值得知道的注意事項。**

- **首次執行成本。** 從 KEGG 冷取得目前所列出的全部約 573 個模組，在 334 ms 的請求間節流（3 req/s）下約需 6–12 分鐘。
  模組 ID 範圍為 `M00001`–`M01063`，但 KEGG 讓此範圍保持稀疏——退役的 ID 不會被重用，因此實際數量低於上界。
  第三階段設定畫面會顯示一個就地進度條，其 ETA 在前 5 個模組完成後，由逐模組實際時間的滾動平均推導而來。
  後續執行會使用快取，且 `Run Enrichment` 按鈕會在數秒內啟用。
- **Group 1 只有兩個選項**（Prokaryotes / Eukaryotes），這在生物學上非常粗略。
  它的存在是為了完整性——例如「任何原核生物」的比較研究——但多數分析會受益於 Level 2（6 個候選）或 Level 3（數十個候選）以取得更細的範圍界定。
- **`min_group_overlap` 是一個研究旋鈕。** 預設 `1`（寬鬆的 ∃-重疊）適合探索性工作。
  對於論文，請考慮測試一個較高的門檻以確保穩健性——一個 Group（例如「Animals」）中數百個生物裡只有一個擁有的模組，即使它通過了預設過濾，對該分析框架而言在生物學上仍是邊緣的。
- **模組 CSV 欄名與路徑模式 CSV 一致。** 兩種模式都匯出相同的標頭：`EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,FDR,HitKeggIDs`。
  （`Expected` 與 `EnrichmentRatio` 的定義見上方逐路徑超幾何步驟：`EnrichmentRatio` 即 fold enrichment = 觀測／期望。）
  在模組模式中，`EntryID` 欄攜帶 `M00001` 式的模組 ID；在路徑模式中則攜帶 `<species_code><pathway_number>` 形式的 ID（例如 `gmx00010`）。

<a id="starting-a-new-analysis-round"></a>
### 開始新一輪分析（Starting a new analysis round）

當你完成某份資料集、想重新開始時，**Start a new analysis** 會清除一切；相對地，步進器的 **Input** 步驟則保留你的設定與快取，讓你能重跑*同一份*資料集。

當你完成一次富集分析、想分析另一份資料集——或從頭重跑整個流程——時，第三階段 **Enrichment Result** 畫面會在 `[Download enrichment results CSV]` 下方獨立一行提供一個 **Start a new analysis** 按鈕。
按下它會開啟一個確認對話框，警告目前的 DAM / 富集分析結果、以及任何尚未下載的圖或 CSV 都將遺失。
按下 **Start over** 後，應用程式會把每個參數重置為其預設值、清除已載入的 MS-DIAL `.txt` / metadata `.csv` 以及記憶體中的 KEGG 資料，並把你帶回第一階段——*而不會*重新執行啟動時的生物清單載入。
（磁碟上的 KEGG 快取會留存，因此事後重新取得相同物種或模組會是快速的快取命中。）

這刻意與階段步進器的 **Input** 步驟有所區別，後者會導航回第一階段，同時*保留*每項設定、已載入的檔案與已取得的快取，讓你能在**同一份**資料集上持續迭代。
用步進器來調整並重跑目前的分析；用 **Start a new analysis** 來捨棄一切並重新開始。
若你之後可能想再用目前的設定，請在重新開始之前透過 Data 分頁的 **[Save settings…]** 按鈕儲存它。

---

<a id="advanced-topics--reference"></a>
## 進階主題與參考（Advanced topics & reference）

其餘各節是你可依需要閱讀的參考材料——基礎的「缺失 vs 零」概念、雙模式輸入、快取、設定檔、問題回報，以及引用文獻。

<a id="missing-values-nan-vs-true-zeros-00"></a>
## 缺失值（`NaN`）與真正的零（`0.0`）（Missing values vs true zeros）

空白儲存格與量測到的零*並不*是同一回事，而 metabolopan 拒絕把它們混為一談。
空白意味著「我們從未量測這個」（`NaN`，內部的*缺失*標記）；`0.0` 則意味著「我們量測了它，而它確實是零」。
本節會明確說明每一種如何被處理，因為「把空白填補為 `0`」這個常見捷徑會默默地使下游每一項統計產生偏誤。

metabolopan 在**一筆缺席的量測**與**一筆確實為零的量測**之間劃出一條明確的界線，而這項區別會被刻意且一致地貫穿每一個下游步驟。
在分析前把缺失儲存格填補為 `0`（一個常見習慣）會默默使統計產生偏誤，所以本節明確說明每個值的意義以及它如何被處理。

**規則，在載入時即固定。** 當解析 MS-DIAL `.txt` 時，一個空白／僅含空白字元／`"null"`／`"NA"`（不分大小寫）或其他無法解析的強度儲存格會變成 `f64::NAN`——*缺失／未量測／無法計算*的內部標記（這裡 `f64::NAN` 即 IEEE 浮點數的「非數值」）。
一個字面上讀作 `0` 的儲存格則解析為真正的 `0.0`。
數值型中繼資料欄位遵循相同的劃分：空白儲存格是*缺席*（`None`），而寫成 `0` 則是真正的零（而且因為它會被當作正規化的除數，會被視為資料輸入錯誤而報錯，而非默默當成缺席）。

**核心行為：`NaN` 被略過、`0.0` 參與運算。** DAM 中每個逐特徵的歸約——群組平均、中位數、變異數、IQR 與相異值計數——都會先*丟棄* `NaN` 值，再對剩下的部分計算。
`0.0` 作為一筆真實觀測，會完整地進入算術。
在同樣三個重複樣本上，這就是差異：

| Group values    | Effective *n* | Mean  | Variance             |
| --------------- | ------------- | ----- | -------------------- |
| `[10, 12, NaN]` | 2             | 11.0  | computed on 2 points |
| `[10, 12, 0]`   | 3             | 7.33  | much larger          |

一個缺失的重複樣本表現得像那個樣本不存在；一個為零的重複樣本則把平均拉低、灌大離散程度，並計入樣本大小。

**這項區別在哪裡浮現，逐步說明：**

| Step | `NaN`（缺失）行為 | `0.0`（真正的零）行為 |
| --- | --- | --- |
| **逐群組前置過濾** | 一個 `NaN` 會降低非 `NaN` 計數；一個完全是 `NaN` 的群組會使該特徵無法檢定，因此它會被略過並計入 `skipped`，而非占用一個「不顯著」的名額。（DAM 要求*每個*群組中**至少有 2 個非 `NaN` 值**。） | 一個 `0.0` 算作存在，且能幫助某群組達到最低數量。（一個全是相同零值的群組則改由「無變異」檢查 `nunique > 1` 與 `IQR > 0` 移除——相同的略過、不同的理由。） |
| **統計檢定** | Student / Welch / Brunner–Munzel 各自計數非 `NaN` 值，並在某群組於 `NaN`-丟棄後少於 2 個時回傳 `NaN` *p* 值。 | 一個 `0.0` 會像任何其他數字一樣流入 *t* 統計量、變異數、標準誤與自由度。 |
| **對數轉換** | `NaN` 會原封不動地通過（轉換會略過它；絕不報錯）。 | 可選的變異數穩定化轉換是 **`arcsinh`（`asinh`），而非 `log10`**——之所以這樣選，是為了讓零是安全的：`asinh(0) = 0`（一個有限、可用的值），而 `log10(0) = −∞`。偏好 `arcsinh` 的刻意理由：無須任何偽計數或夾限即可承受零。 |
| **倍數變化** | 一個 `NaN` 無法驅動 `±∞` 倍數變化，因為它一開始就被排除於平均之外——`NaN` 倍數變化是保留給「該值確實無法被計算」的情形。 | 只有真正的 `0.0` 才能把某組的平均（或中位數）壓到剛好為 0，從而使 `log2(FC) = ±∞`；這些特徵會被停靠在火山圖 X 軸的邊緣（加上小幅 jitter），絕不會被默默丟棄。 |
| **FDR 校正** | Benjamini–Hochberg / Benjamini–Yekutieli 會完全略過 `NaN` *p* 值：它們**不會**占用被校正的 *m* 個檢定之一，且 `NaN` 會原封不動地通過到輸出。 | 一個有限的 *p* 值——包括由資料中真正的零所產生的——則會正常被校正。 |
| **趨勢分類** | 一個其校正後 *p* 值或 `log2(FC)` 為 `NaN` 的特徵會被分類為 `NotSignificant`；它永遠不會被判為 Up 或 Down。 | 一個有限、顯著的結果會照常被分類為 Up 或 Down。 |
| **`NaN` 與 `±∞` 保持區別** | `NaN` 意味著「無法被計算」（一個 *n* < 2 的群組；Brunner–Munzel 下完美分層的群組）。`NaN` 點會從圖中丟棄，但仍會列於 CSV 中。 | `±∞`（即正負無限大）是一個真實、*有序*的結果——一個 underflow 到剛好為 0 的 *q* 值在 `−log10` 軸上會變成 `+∞`（*超出刻度但有序*），而一個零平均群組會給出 `±∞` 倍數變化。`+∞` 點會被停靠在圖的邊緣。 |

> **注意：** `f64::NAN` 是*缺失*標記，`f64::INFINITY` / `±∞`（即正負無限大）是一個*超出刻度但有序*的真實結果，而 `0.0` 是*一個真正的零*——這三個不同的狀態，軟體絕不會把任何一個塌縮成另一個。

**CSV 編碼。** 匯出時每個狀態都會被寫成不同形式，使檔案可往返：

| Value | Written as |
| --- | --- |
| `NaN` | empty（`""`）——若檔案被重新讀回，會還原為「缺失」 |
| `+∞` | `inf` |
| `−∞` | `-inf` |
| `0.0` | `0` |

**一個刻意的例外。** PQN 正規化把每特徵參考商為 `0` 的情形與 `NaN` 同等對待：兩者都被排除為*不可用*，因為零商對 PQN 的「商之中位數」因子不攜帶任何資訊。
這是唯一一處兩者被刻意合併的地方。

**結論。** 把確實缺失的量測留為*空白*，並把量測到的零寫成 `0`；軟體會從輸入到匯出都讓它們保持區別。
若你事先把缺失儲存格填補為 `0`，就會灌大樣本數、把群組平均拉向零、扭曲變異數與倍數變化，並使差異累積的判定偏誤——所以請讓 metabolopan 以 `NaN` 承載缺失值，由它替你完成這些記帳。
<a id="dual-mode-positive--negative-ionization-input"></a>
## 雙模式（正 + 負離子）輸入（Dual-mode input）

若你把同一批樣本同時跑過正離子化與負離子化，你會有兩個描述同一個實驗的 MS-DIAL `.txt` 檔。
雙模式會一次載入兩者，並以一條刻意保守的聯集規則融合它們的富集訊號——只有在沒有任何模式持反對意見時，某化合物才算作「up」。
中繼資料中的 `biosample` 欄，正是用來告訴工具 `CTR_positive_01` 與 `CTR_negative_01` 其實是同一個生物重複。

代謝體學實驗常把同一批生物樣本同時跑正離子化與負離子化兩種模式，每個研究因此產生兩個 MS-DIAL `.txt` 匯出檔。
本應用程式支援一次載入兩個檔案，並以一條保守的聯集規則合併它們的富集訊號。

<a id="when-to-use-dual-mode"></a>
### 何時使用雙模式（When to use dual-mode）

只要你對同一批生物樣本同時擁有 POS 與 NEG 兩個 `.txt`、並想要一份能反映「任一離子化所提供證據」的單一富集結果，就使用雙模式。
單模式（一個 `.txt`）仍是預設。

<a id="preparing-inputs"></a>
### 準備輸入（Preparing inputs）

1. **兩個 `.txt` 檔。** 每個離子模式一個。
   `Adduct type` 欄同時驅動兩件事：插槽 1 模式單選按鈕的自動填入（見下方 *第一階段 UI*），以及當使用者手動覆寫為相反極性時的不一致提示（以 `+` 結尾的加成物推論為 Positive、以 `-` 結尾推論為 Negative）。
2. **一個含 `biosample` 欄的中繼資料 CSV**（例如標頭 `sample,biosample,group`；欄位順序不拘）。
   每一列把一個逐模式的樣本名稱（例如 `CTR_positive_01`、`CTR_negative_01`）對應到其**生物樣本標籤**（兩個模式皆為同一個 `CTR-01`）與群組。
   biosample 欄讓工具能辨識兩個名稱不同的樣本其實是同一個生物重複。

以沒有 `biosample` 欄的 CSV 進行的雙模式執行，會在第一階段被一個明確的錯誤擋下——請新增 `biosample` 欄、或移除第二個 `.txt` 以繼續。

> **單模式並不需要 `biosample` 欄。** 它只在載入第二個 `.txt` 時才為必需。只有一個 `.txt` 時，單純的 `sample,group` 形式就足夠
> （若存在 `biosample` 欄，會以名稱辨識並排除於第二階段中繼資料正規化的單選按鈕之外——它不會被當作數值中繼資料欄提供）。

<a id="unbalanced-or-missing-mode-samples"></a>
### 不平衡或缺少某模式的樣本（Unbalanced or missing-mode samples）

下方的計算範例是完全平衡的（每個生物樣本都在兩個模式中跑過），但真實研究有時只在單一極性下採集某個生物樣本。
`biosample` 欄是配對兩個模式的依據，因此第一階段會在允許 `Continue to DAM` **之前**強制執行三項雙模式完整性檢查。
每一項都會浮現一個明確的錯誤：

1. **每個群組在*每個*模式中都需要 ≥ 2 個樣本。** 若某群組在 POS 有足夠的重複、但在 NEG 掉到 2 以下（例如因為數個生物樣本缺少 NEG 採集），第一階段會擋下並顯示 `Group 'X' has N sample(s) in POS but M in NEG — both modes need ≥ 2.`。只要每個群組仍各自跨過「每模式 2 個」的門檻，少數缺某模式的樣本是可容忍的；真正把關的只有「逐群組、逐模式」的計數。
2. **生物樣本在同一模式內必須唯一。** 兩列把同一個生物樣本標籤對應到同一個模式，會觸發 `Biosample 'B' appears in N POS rows — must be unique per mode.`。
3. **生物樣本在跨模式間必須維持同一群組。** 若 `CTR-01` 在 POS 是 `control`、在 NEG 卻是 `treatment`，第一階段會擋下並顯示 `Biosample 'B' is in group 'X' in POS but 'Y' in NEG.`。

**這三則訊息中的 `POS` / `NEG` 標籤跟隨各插槽，而非固定順序。** 每個標籤就是該插槽實際被設定的模式，依插槽順序讀取（先插槽 1、再插槽 2）。
上述範例假設第一階段自動填入的常見配置：插槽 1 = Positive、插槽 2 = Negative；若把 Negative 放在插槽 1，同樣的錯誤會以對調的模式讀出（例如 `… N sample(s) in NEG but M in POS …`）。

**通過關卡的缺某模式樣本所造成的影響。** 兩個模式各自在自己的樣本欄上獨立執行 DAM——某個在 NEG 缺席的生物樣本，單純就不會在 NEG 執行中被迭代，因此該模式對其群組的重複較少、檢定力也相應較低；但這不會使該次執行失效。
在第三階段，聯集是在**化合物**層級建立（依下方「僅衝突即嚴格」規則）、而非樣本層級，因此某生物樣本缺少一個模式，只會讓該模式對受影響的化合物貢獻 `Absent`——整合後的 K 不受影響。

**建議。** 為了得到最乾淨的雙模式結果，請在兩種極性下都採集每個生物樣本。
若某些樣本確實只屬於單一極性，請將它們只保留在其存在的那個模式中（只要每個群組仍各自有「每模式 ≥ 2 個」），或捨棄不平衡的那一側。
在 `.txt` 中出現、但不在中繼資料 CSV 中的樣本，會被標記為 `Unassigned`，並在第一階段 → 第二階段的邊界自動捨棄（見上方的群組對應說明），這是另一種排除不需要欄位的方式。

<a id="stage-1-ui"></a>
### 第一階段 UI（Stage 1 UI）

插槽 #1（永遠可見）與插槽 #2（由 `+ Add second ionization mode` 按鈕展開）各有一個檔案選擇器、一個模式單選按鈕（Positive / Negative）、以及一個逐插槽摘要。
插槽 1 的模式單選按鈕會在每次重新載入檔案及重新挑選時，依 `infer_polarity(&table)` 自動填入：`≥ 95%` 為正極性後綴的 Adduct 欄會設為 Positive、`≥ 95%` 為負極性後綴設為 Negative、含糊的混合則讓單選按鈕保持未設定（既有的灰色「Could not auto-detect…」提示仍然適用）。
當插槽 1 的模式被設定後，插槽 2 的單選按鈕會在三種觸發時機自動填入為**相反**值：(1) 透過 `+ Add second ionization mode` 按鈕展開插槽 2、(2) 載入插槽 2 的 `.txt`、(3) 插槽 1 的模式變更（手動點擊或重新挑選重新推論）——情況 (3) 在新的插槽 1 值與插槽 2 已顯示的值衝突時，也會翻轉插槽 2。
使用者仍可手動點擊任一單選按鈕來覆寫。
插槽 2 的單選按鈕仍會停用已被插槽 #1 選走的那個選項（工具提示會說明原因）。
加成物不一致提示（「黃色：Adduct column says X but you selected Y」）仍會在手動覆寫與自動偵測相左時觸發；兩種提示都不會擋下 `Continue to DAM`。

<a id="stage-2-shared-setup-per-mode-dam"></a>
### 第二階段（共用設定、逐模式 DAM）（Stage 2）

第二階段使用單一設定畫面——一種正規化方法、一組比較（分子組 vs 分母組）、一種 DAM 方法、一種 FDR 方法——並套用於**兩個**模式。
在執行器內部，兩個 tokio worker 會並行對每個模式執行 `run_dam`；執行畫面會顯示兩條堆疊的進度條。
若任一模式失敗，錯誤訊息會指出是哪個模式（`POS: ...` 或 `NEG: ...`）。

火山圖畫面會在繪圖區上方繪製一條 `POS | NEG` 分頁列。每個分頁各自快取自己的紋理；變更任一門檻滑桿都會使兩者失效。
PNG 匯出使用逐模式的預設檔名（`volcano-pos.png` / `volcano-neg.png`）。
DAM CSV 匯出會在開頭發出一行 `# Mode: dual (POS+NEG)` 註解，並在每一列前面附加一個 `Mode` 欄，列順序為先 POS、後 NEG。

<a id="stage-3--dual-mode-n-and-k-math"></a>
### 第三階段 — 雙模式 N 與 K 的運算（Stage 3 — dual-mode N and K math）

第三階段在「僅衝突即嚴格」聯集規則下，從**兩個**模式的 DAM 特徵建立母體 N 與前景 K（保守的選擇：方向相反的訊號會排除某化合物）。

**PubChem 與 KEGG `/conv` 呼叫只在聯集後的 InChIKey 集合上執行一次**，使網路成本在雙模式下不會加倍。

**N（母體）** = 經由 PubChem → KEGG conv 鏈，從任一模式任一特徵可達之每個 cpd 的聯集。

**逐模式趨勢彙總。** 對每個 cpd `c`，分別從每個模式收集逐特徵趨勢並加以彙總。五種可能的逐模式判定為：

| Trend | 此 cpd 在此模式中的意義 |
| --- | --- |
| `Up` | 此模式中有任何貢獻特徵被標記為 Up，且無任何 Down |
| `Down` | 對稱情形（有任何 Down，且無任何 Up） |
| `NS` | 只有不顯著的特徵 |
| `Conflict` | 同一模式中同時有 Up 與 Down 特徵（same-InChIKey-different-trends 的邊界情況） |
| `Absent` | 此 cpd 從此模式完全不可達 |

**「僅衝突即嚴格」規則下的 K（前景）。** 對於作用方向 `Up`：某 cpd 進入 K 的條件為——至少一個模式說 Up，且沒有任何模式說 Down，且沒有任何模式處於 Conflict。
`Down` 為對稱情形。
`Both` 要求至少有一個 Up 或 Down 訊號，且無任何 Conflict，且非（一個模式 Up 而另一個模式 Down）。

**單模式套用相同的衝突規則。** 單模式執行是此規則的退化單模式情形：某化合物同時被該單一模式內的一個 Up 特徵與一個 Down 特徵到達——兩個不同的 InChIKey 對應到**同一個** KEGG 化合物，一個 Up + 一個 Down——會彙總為 `Conflict` 並被**排除**於 K 之外，與雙模式採取相同的保守選擇。
（在此之前，單模式會把這類含糊化合物保留在 K 中。）因衝突被排除的計數會出現在第三階段的 INFO 日誌中。
對於任何沒有這種模式內衝突的資料集，單模式的 K 維持不變。

底部面板的 **Data** 分頁會把雙模式分割呈現為母體 / 前景來源漏斗的一部分：

```
Universe — all tested features (measurable metabolome)
  … → N KEGG cpds  (POS-only: a; NEG-only: b; in both: c)
Foreground — significant features (active direction)
  … → K KEGG cpds  (sig POS-only: d; sig NEG-only: e; agree both: f; excluded by conflict: g)
```

當某個模式貢獻了每一個 K cpd 時，會出現一行黃色的 `K source: POS only (NEG had 0 sig features in the active direction)`。
富集 CSV 會在開頭發出一行 `# Mode: dual (POS+NEG)` 註解；逐列的 CSV 形狀不變（ORA 運算與模式無關）。

<a id="worked-example"></a>
### 計算範例（Worked example）

使用 `data/double-mode/` 的測試固定資料（8 Treatment + 8 Control + 3 QC 個生物樣本，每個都在兩個模式中採集 = 跨 19 個生物樣本的 38 個樣本欄；中繼資料另帶一個數值 `mass` 欄）：

1. 第一階段：把 `data-positive.txt` 載入插槽 #1（Mode: Positive）、把 `data-negative.txt` 載入插槽 #2（Mode: Negative），並載入 `metadata.csv`。
   點擊 `Continue to DAM`。
2. 第二階段：挑選 `Treatment` vs `Control`（第三個 `QC` 群組不選），正規化與 FDR 維持預設。
   執行畫面會顯示兩條進度條；每個模式約需 6–60 秒，視特徵數而定。
3. 第二階段門檻：在 POS 與 NEG 分頁間切換以檢視各自的火山圖；下載分頁式 PNG 或統一的 CSV。
   點擊 `Continue to Enrichment`。
4. 第三階段設定：挑選一個 KEGG 物種（路徑模式）或 Level + 生物群組（模組模式）；行內進度列會串流 KEGG 擷取進度。
   完成後，點擊 `Run Enrichment`。
5. 第三階段結果：結果面板會顯示分解區塊；因衝突被排除的 cpd ID 會以 INFO 等級出現在日誌中。
   若想在點圖上少一些或多一些列，可行內調整 Top N，再點擊 `Re-draw dot plot`。
   點圖會保留 Top N 個*最顯著*的條目，並依*富集倍數*堆疊（最大者在最上方）；畫布高度會在每次重繪時重新貼合所顯示的列數（見上方「點圖的選取 vs 排序」）。

<a id="caches-and-provenance"></a>
## 快取與來源（Caches and provenance）

為避免每次工作階段都重新下載相同的 KEGG / PubChem 資料，本應用程式會保留本地快取檔且永不使其過期——不論年代多久，已快取的條目都會被回傳，何時刷新由你決定。
**Data** 分頁會中立地顯示每個快取的擷取日期，讓你自行判斷新鮮度。

**磁碟上的檔案**（位於 KEGG 快取目錄中）：

- `inchikey.json` — PubChem InChIKey → CID 結果。
- `cid_to_cpd.json` — KEGG CID → 化合物結果。
- `modules.json` — 已擷取的 KEGG 模組條目。
- `organisms.json` — KEGG 生物名冊（啟動時載入一次）。
- `.inchikey.lock` / `.cid_to_cpd.lock` — 短命的寫入鎖（以點為前綴／隱藏）。
- `.modules.lock` — 長時間運行的模組擷取建議鎖。

第三階段的快取（`inchikey.json`、`cid_to_cpd.json`、`modules.json`）儲存**逐條目**的 `fetched_at` 時間戳，有別於第一階段的物種快取（檔案層級時間戳）。
逐條目的粒度是刻意的：這些快取會在數週或數月、橫跨許多工作階段中增量成長，而檔案層級時間戳要嘛對年代撒謊、要嘛強迫頻繁的完整刷新。
第三階段結果畫面會把這呈現為一個時間跨度（`PubChem CIDs fetched date: 2026-03-01 -> 2026-05-22 (<n> entries used)`）；模組模式還會額外顯示模組快取在該次執行所用之**保留**模組上的時間跨度，而非整個快取。

> **給開發者：** 每個逐條目的 `fetched_at` 都是一個 `DateTime<Utc>`（一個 UTC 時間戳）。
> 快取鎖機制：
> - **PubChem `.inchikey.lock` + KEGG `.cid_to_cpd.lock`** — 短命，只在快取寫入期間持有。等待 30 s、以 100 ms 輪詢。（兩個檔案皆以點為前綴／隱藏。）
> - **KEGG `.modules.lock`** — 長時間運行的建議鎖，在整個約 6–12 分鐘的模組擷取期間持有。鎖檔內嵌持有者的 PID 與一個至多每 30 s 改寫一次的心跳 `last_seen_at` 時間戳。並行的應用程式實例會看見這個存活的鎖，並等待至多 30 min（5 s 輪詢）直到它清除。若心跳超過 90 s 未更新，該鎖會被視為孤立（持有者已崩潰）並被覆寫。這可防止兩個應用程式實例在模組擷取迴圈中競速、進而一同觸發 KEGG 的 403 速率限制。
> - **啟動清理。** 每次應用程式啟動時，快取目錄的鎖檔（`.inchikey.lock`、`.cid_to_cpd.lock`、`.modules.lock`）都會被無條件移除，使一次崩潰絕不會永久阻擋未來的寫入。

快取新鮮度——**沒有過期門檻**。
所有 KEGG 快取都不會過期：不論年代多久，已快取的條目都會被回傳，且應用程式絕不會自行默默重新擷取。
取而代之的是底部面板 **Data** 分頁的 `Cache data` 區塊（位於 Enrichment Analysis + Enrichment Result 畫面）會中立地呈現擷取時間，把刷新的決定留給你：

- 逐物種路徑快取：顯示 `KEGG pathways (<code>): <ts>`（在設定與結果畫面皆有）；經由 `Refresh KEGG pathway cache` 按鈕重新擷取。
- 模組條目：顯示一個 `KEGG modules fetched date: <oldest> -> <newest>` 跨度；warm-fetch 的決定取決於快取鍵的成員資格。
  經由 `Refresh KEGG module cache` 按鈕重新擷取。
- 在 Enrichment **Result** 畫面上，目錄刷新按鈕（模組 / 路徑）會導回設定畫面，在那裡執行重新擷取（其進度列位於該處）；PubChem / KEGG-conv 的刷新則透過一個確認對話框就地執行。
- 生物名冊（`organisms.json`）：啟動時載入一次（快取優先：不論年代多久，磁碟上的副本永遠勝出），可在應用程式內透過 Data 分頁 `Cache data` 區塊中的 `Refresh KEGG organism list` 按鈕刷新。
  該按鈕會就地重新擷取 `/list/organism` 而無需重新啟動；或者，從快取目錄刪除 `organisms.json` 並重新啟動以強制冷擷取。
  （`Refresh KEGG pathway cache` 按鈕是分開的——它只重新擷取所選物種的「路徑→化合物」對應，而非生物名冊。）

<a id="saving-and-loading-session-settings-reproducibility"></a>
## 儲存與載入工作階段設定（再現性）（Saving and loading session settings）

你通常永遠不需要手動碰這個檔案——當你點擊 **[Save settings…]** 時，應用程式會替你寫入它，並在 **[Load settings…]** 時把它讀回來。
它的存在是為了讓一次執行*可再現*：把同一份快照加上同一份輸入交給合作者（或未來的你），分析結果就會逐位元相同。
此處之所以記載格式，只是給那些想要將它腳本化、或想檢視擷取了什麼的人參考。

Data 分頁中的兩個按鈕——**[Save settings…]** 與 **[Load settings…]**——讓你把每一項第一至第三階段的參數快照到一個 JSON 檔，並在日後重新套用。
其用意在於再現性：若你（或合作者）以同一份快照與同一份輸入重新執行，分析結果會逐位元相同。

<a id="whats-in-the-file"></a>
### 檔案內容（What's in the file）

一份美化排版（pretty-printed）的 JSON，包含：

- `schema_version`（目前為 `1`——磁碟上的 schema 基準）、`app_version`、`saved_at`（UTC）、一個初始為 `""` 的 `user_note` 欄——你可以用任何文字編輯器打開檔案並填入它。
- `input_files` — 對於你在儲存時已載入的每個 MS-DIAL `.txt` 與中繼資料 `.csv`：該檔的 basename + 其 SHA-256。
  **僅雜湊——絕不包含你的原始資料。** 這讓未來的 Load 能偵測到你的輸入是否已偏離當初製作快照所依據的版本。
- `settings` — 從第一階段到第三階段的每一項參數（分析模式、物種 / 生物群組、比較群組、DAM 方法、正規化、FDR 方法、門檻、匯出尺寸、富集方向 / FDR / top-N）。

<a id="the-full-file-field-by-field"></a>
### 完整檔案，逐欄位（The full file, field by field）

一個完整的範例（分析進行到一半時取的單模式快照）。外層欄位已於上方說明；下表記載 `settings` 之下的每一個鍵。
此範例顯示了被填入的選用欄位——`null` 是它們的預設（見下表）。

```json
{
  "schema_version": 1,
  "app_version": "1.2.3",
  "saved_at": "2026-06-04T09:15:22Z",
  "user_note": "",
  "input_files": [
    { "role": "positive", "name": "MS-DIAL-output-example.txt", "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" },
    { "role": "metadata", "name": "metadata-example.csv",       "sha256": "9f86d081884c7d659a2feaa0c55ad015a3bf4f1b2b0b822cd15d6c15b0f00a08" }
  ],
  "settings": {
    "analysis_mode": "Pathway",
    "kegg_species": "hsa",
    "organism_group_level": null,
    "organism_group": null,
    "min_group_overlap": 1,
    "numerator": "Treatment",
    "denominator": "Control",
    "dam_method": "Student",
    "drop_unknown": true,
    "dedup_enabled": true,
    "normalization": "None",
    "metadata_column": null,
    "pqn_reference": "AllSamples",
    "pqn_reference_group": null,
    "log_transform": true,
    "dam_fdr_method": "BenjaminiHochberg",
    "fc_threshold": 2.0,
    "fdr_threshold": 0.05,
    "delta_threshold": 0.33,
    "stage2_export_width_in": 3.5,
    "stage2_export_height_in": 2.2,
    "stage2_export_dpi": 300,
    "direction": "Both",
    "top_n": 20,
    "enrichment_fdr_threshold": 0.05,
    "min_hit_count": 1,
    "min_entry_size": 1,
    "enrichment_fdr_method": "BenjaminiYekutieli",
    "stage3_export_width_in": 3.5,
    "stage3_export_height_in": 7.0,
    "stage3_export_dpi": 300
  }
}
```

外層：`schema_version` 必須為 `1`（其他值在 Load 時會被拒絕）；`app_version` / `saved_at` 為提示性資訊；`user_note` 是你可手動編輯的自由文字；每個 `input_files` 條目是 `role`（`positive` / `negative` / `metadata`）+ 檔案 basename + SHA-256（僅雜湊——絕不含原始資料）。

> **注意：** 少數鍵使用*物件變體*（object-variant）語法——值不是裸字串，而是一個攜帶資料的小物件，例如中繼資料正規化的 `{"Metadata":{"column":"<name>"}}`、或逐群組 PQN 參考的 `{"Group":"<name>"}`。外層鍵（`Metadata`、`Group`）為變體命名；內層物件持有其參數。

`settings` 之下的每一個鍵。**UI control** 欄把每個鍵對應到設定它的畫面／控制項：

| Key | JSON type / allowed values | Default | UI control | 意義與限制 |
| --- | --- | --- | --- | --- |
| `analysis_mode` | `"Pathway"` \| `"Module"` | `"Pathway"` | **Analysis Mode** 單選按鈕（第三階段設定） | 第三階段 ORA 目錄：逐物種路徑 vs 逐 Group 模組。 |
| `kegg_species` | string \| `null` | `null` | KEGG 物種選擇器（第三階段設定，路徑模式） | 路徑模式用的 KEGG 生物代碼（例如 `"hsa"`）。 |
| `organism_group_level` | `1`–`3` \| `null` | `null` | Level 單選按鈕（第三階段設定，模組模式） | KEGG 生物階層層級（模組模式）。 |
| `organism_group` | string \| `null` | `null` | Group 下拉選單（第三階段設定，模組模式） | 所選的生物 Group 名稱（模組模式）。 |
| `min_group_overlap` | integer ≥ `1` | `1` | **Minimum group overlap** 控制項（第三階段設定，模組模式） | 模組模式：只在某模組與所選 Group 共享 ≥ 此數的生物時才納入。透過第三階段設定畫面的 **Minimum group overlap** 控制項設定（範圍 `1`–`min(Group size, 20)`）；也記錄於匯出 CSV 的 `# MinGroupOverlap:` 行。 |
| `numerator` | string \| `null` | `null` | Numerator group ComboBox（第二階段設定） | DAM 分子組。若不存在於目前的中繼資料中，Load 時重設為 `null`。 |
| `denominator` | string \| `null` | `null` | Denominator group ComboBox（第二階段設定） | DAM 分母組；必須與 `numerator` 不同（於 Start DAM 時檢查）。若不存在於目前的中繼資料中，Load 時重設為 `null`。 |
| `dam_method` | `"Student"` \| `"Welch"` \| `"BrunnerMunzel"` | `"Student"` | **DAM method** 單選按鈕（第二階段設定） | DAM 統計檢定。 |
| `drop_unknown` | `true` \| `false` | `true` | **Drop Unknown** 切換（第二階段設定） | 在檢定前捨棄 InChIKey 為 null 的特徵。 |
| `dedup_enabled` | `true` \| `false` | `true` | **Dedup** 切換（第二階段設定） | 以 InChIKey 去除重複特徵（級聯）。 |
| `normalization` | `"None"` \| `"Sum"` \| `"Median"` \| `"Quantile"` \| `{"Metadata":{"column":"<name>"}}` \| `{"Pqn":{"reference":<pqn_reference>}}` | `"None"` | Normalization 單選按鈕（第二階段設定） | 樣本軸正規化。`Metadata` 與 `Pqn` 是攜帶資料的物件變體。 |
| `metadata_column` | string \| `null` | `null` | Metadata-column ComboBox（第二階段設定，Metadata 正規化） | `Metadata` 正規化所用的欄。若在目前資料中不是數值中繼資料欄，Load 時重設為 `null`。 |
| `pqn_reference` | `"AllSamples"` \| `{"Group":"<name>"}` | `"AllSamples"` | PQN reference 單選按鈕（第二階段設定，PQN 正規化） | PQN 參考光譜（只在 `normalization` 為 `Pqn` 時有意義）。 |
| `pqn_reference_group` | string \| `null` | `null` | PQN reference-group ComboBox（第二階段設定） | 當 `pqn_reference` 為 `{"Group":…}` 時的 Group 名稱。若不存在於目前的中繼資料中，Load 時重設為 `null`。 |
| `log_transform` | `true` \| `false` | `true` | **Log transformation** 切換（第二階段設定） | 在 Welch/Student 之前套用 arcsinh（BM 會忽略它）。當手動編輯的 v1 檔缺少此鍵時，預設為 `true`。 |
| `dam_fdr_method` | `"BenjaminiHochberg"` \| `"BenjaminiYekutieli"` | `"BenjaminiHochberg"` | 第二階段 FDR 單選按鈕（第二階段設定） | 第二階段 FDR。`"NoCorrection"` 在 **Load 時會被強制轉為 BH**（第二階段絕不暴露 None）。 |
| `fc_threshold` | `1.0`–`1024.0` | `2.0` | Fold-change 門檻（第二階段結果） | 火山圖／CSV 的倍數變化截斷（使用 `\|log2(FC)\| > log2(value)`）。 |
| `fdr_threshold` | `0.0001`–`1.0` | `0.05` | FDR 門檻（第二階段結果） | 火山圖／CSV 的 q 值截斷。 |
| `delta_threshold` | `0.0`–`1.0` | `0.33` | Cliff's δ 門檻（第二階段結果，僅 BM） | Cliff's δ 截斷（僅 Brunner–Munzel；Welch/Student 會忽略）。 |
| `stage2_export_width_in` | `1.0`–`40.0` | `3.5` | **Width (in)** 欄位（第二階段結果） | 火山圖 PNG 寬度（英吋）。 |
| `stage2_export_height_in` | `1.0`–`40.0` | `2.2` | **Height (in)** 欄位（第二階段結果） | 火山圖 PNG 高度（英吋）。 |
| `stage2_export_dpi` | `72`–`1200` | `300` | **DPI** 欄位（第二階段結果） | 火山圖 PNG 解析度。 |
| `direction` | `"Up"` \| `"Down"` \| `"Both"` | `"Both"` | **Include DAM features with direction** 單選按鈕（第三階段設定） | 哪些 DAM 特徵組成 ORA 前景（UI：Up only / Down only / Both）。 |
| `top_n` | `1`–`100` | `20` | **Top N pathways** 輸入（第三階段結果） | 點圖上繪製的最大條目數。 |
| `enrichment_fdr_threshold` | `0.0001`–`1.0` | `0.05` | **Enrichment FDR threshold**（第三階段結果） | ORA 顯示的顯著性截斷。 |
| `min_hit_count` | `1`–`10` | `1` | **Minimum hit count**（第三階段結果） | FDR 後的顯示過濾：隱藏命中數較少的條目。 |
| `min_entry_size` | `1`–`20` | `1` | **Minimum number of compounds detected in a pathway/module**（第三階段設定） | FDR 前的條目過濾：捨棄母體化合物數少於此值的條目。當手動編輯的 v1 檔缺少此鍵時，預設為 `1`。 |
| `enrichment_fdr_method` | `"BenjaminiHochberg"` \| `"BenjaminiYekutieli"` \| `"NoCorrection"` | `"BenjaminiYekutieli"` | **FDR correction** 單選按鈕（第三階段設定） | 第三階段 FDR（預設 BY——ORA 條目共享化合物）。此處允許 `"NoCorrection"`（與第二階段不同）。 |
| `stage3_export_width_in` | `1.0`–`40.0` | `3.5` | **Width (in)** 欄位（第三階段結果） | 點圖 PNG 寬度（英吋）。 |
| `stage3_export_height_in` | `1.0`–`40.0` | `7.0` | **Height (in)** 欄位（第三階段結果） | 點圖 PNG 高度（英吋）；除非被覆寫，否則自動貼合列數（見第三階段路徑模式之 *11. Exporting the dot plot*）。 |
| `stage3_export_dpi` | `72`–`1200` | `300` | **DPI** 欄位（第三階段結果） | 點圖 PNG 解析度。 |

**這些範圍是應用程式內的控制項上限，而非檔案的硬性限制。** 手動編輯成超出所列範圍的值會照寫載入，只在你下次於應用程式中碰觸該控制項時才會被夾限；匯出尺寸還會額外被夾限，使 `round(inches × DPI)` 在算繪時每軸維持在 `64–20000` px 內。拼錯或多餘的鍵在 Load 時會被拒絕（檔案必須恰好包含這些鍵），任何非 `1` 的 `schema_version` 亦然。上述四個依賴輸入的欄位是唯一會在 Load 時被重設的欄位。

<a id="when-is-each-button-available"></a>
### 各按鈕何時可用（When is each button available）

- **Save settings…** 在啟動畫面之後的每個畫面上都啟用，無論是否已載入輸入。
  從空白的第一階段儲存，會把你偏好的預設擷取為下次可用的預設組合。
- **Load settings…** **只在第一階段**啟用。
  在其他階段該按鈕會變灰；停留滑鼠其上會顯示「Loading settings is only available on the Stage 1 input screen.」這是刻意的——在分析進行到一半時套用快照，會使螢幕上的結果與新參數不同步，因此工作流程要求你從輸入重新執行。

<a id="loading-workflow"></a>
### 載入流程（Loading workflow）

1. 在第一階段點擊 **[Load settings…]**。
   作業系統的檔案選擇器會開啟。
2. 挑選一個已儲存的 `.json`。
   一個確認對話框會向你顯示其中的內容：
   - 儲存時間戳（以你的當地時間顯示）、快照的應用程式版本、使用者備註（若有）。
   - 設定的單行摘要（分析模式、DAM 方法 + FDR、正規化、富集方向 + FDR + top-N）。
   - **雜湊不符**——若你目前載入的任一輸入檔的 SHA-256 與快照不同，會列在此處。
     若你選擇繼續，設定仍會套用，但你會被警告輸入已偏移。
   - **欄位重設**——若快照指名了一個分子 / 分母組、一個中繼資料欄、或一個 PQN 參考群組，而它不存在於你目前載入的中繼資料中，這些欄位會被列出並在套用時重設為 `None`。
     你會需要在第二階段設定中重新挑選它們。
     （此區段只在你於 Load 時已載入中繼資料時出現；若你在上傳中繼資料前就 Load，則安全網改為第二階段設定的關卡——見下一段。）
3. 點擊 **Apply settings** 以覆寫你目前的設定，或 **Cancel** 以丟棄。

<a id="what-if-i-load-settings-before-uploading-metadata"></a>
### 若我在上傳中繼資料前就載入設定，會怎樣？（What if I load settings before uploading metadata?）

快照的 `numerator` / `denominator` 會被原樣寫入設定中（Load 時不進行驗證，因為沒有中繼資料可供比對）。
當你稍後上傳中繼資料並前進到第二階段設定時，關卡會檢查群組成員資格：若保留的值未出現在新中繼資料的群組中，「Start DAM」按鈕會被算繪為灰色，並附帶一行行內警告（`⚠ Numerator/denominator group not present in the loaded metadata.`）以及作為停留提示的相同文字。
從 ComboBox 下拉選單重新挑選一個有效的群組，警告即會清除。

<a id="hand-editing-the-json"></a>
### 手動編輯 JSON（Hand-editing the JSON）

該檔案是純 UTF-8 JSON，美化排版過。
你可以：

- 在 `user_note` 欄加入備註。
- 微調單一門檻而無需從應用程式重新儲存。
- 移除 `input_files` 區塊以分享一份「僅設定」的快照（Load 能處理空的 `input_files` 陣列——雜湊檢查會被略過）。

把 `schema_version` 手動編輯成非 `1` 的數字、或破壞 JSON 語法，會在 Load 時浮現一則明確的錯誤 toast（例如
*"This settings file uses schema version 2; this app expects version 1."* 或 *"Settings file is not valid JSON (line 7 column 15) …"*）。
任何攜帶非 `1` 之 `schema_version` 的快照都會被拒絕——請從你目前的設定重新儲存，以產生一份 v1 快照。

<a id="reporting-bugs"></a>
## 回報問題（Reporting bugs）

若有什麼看起來不對勁——一個錯誤、一次卡死、兜不攏的結果——取得協助最快的方式是點擊日誌窗格中的 **[Download bug report…]**，並把產生的 zip 附到 GitHub issue 或電子郵件。
這個套組在設計上即受隱私邊界約束：它攜帶日誌與設定，絕不含你的原始資料，並從任何路徑中清除你的家目錄。

若有什麼看起來不對勁——一個非預期的錯誤、一個卡死的階段、與預期不符的結果——取得協助最簡單的方式是點擊日誌窗格中的 **[Download bug report…]**（位於視窗底部、**Clear** 按鈕旁）。
一個確認對話框會列出產生的 zip 將包含的檔案，接著一個儲存檔案對話框讓你挑選要放在哪裡。

該 zip 恰好包含八個檔案：

- `README.txt` — 說明該套組及其隱私邊界。
- `version.txt` — 應用程式建置資訊（套件版本、git SHA、rustc、target）。
- `RUST_LOG.txt` — 僅 `RUST_LOG` 指示詞的值，單獨一行。
- `KEGG_CACHE_DIR.txt` — 僅 `KEGG_CACHE_DIR` 環境變數的值（或 `<unset>`）。
  這兩個是逐變數的檔案（檔名 = 變數名），因此沒有人會把這個套組誤認為完整的環境傾印——只有這兩個具名變數會被納入。
- `logs.txt` — 本工作階段的每一則 INFO / WARN / ERROR 事件（HTTP 與其他低階依賴的雜訊會被過濾掉，使檔案保持可讀）。
- `app_state.txt` — 你所在的階段與目前的設定（分析模式、物種／群組、比較群組、FDR 方法、門檻等）。
- `input_summary.txt` — 你已載入的 MS-DIAL 檔案與中繼資料 CSV 的路徑與計數（僅路徑——無儲存格值）。
- `cache_summary.txt` — KEGG / PubChem 快取檔的大小與新鮮度時間戳（無快取內容）。

**隱私：**

- 此套組絕不包含你的原始 MS-DIAL `.txt` 輸入、你的中繼資料 CSV、或任何先前的 CSV/PNG 匯出。
- 套組內的絕對路徑會把你的家目錄替換為 `~`（例如
  `/Users/alice/Projects/study/POS.txt` 變成 `~/Projects/study/POS.txt`），使該套組在公開分享時（GitHub issue、電子郵件）不會洩漏你的帳號／使用者名稱。
- 只有 `RUST_LOG` 與 `KEGG_CACHE_DIR` 環境變數會被呈現——絕不暴露完整的程序環境。

你可以安心地把這個 zip 附到 GitHub issue 或以電子郵件寄出，無需擔心洩漏你的實驗資料或機器身分。

逐工作階段的日誌檔也會保留在磁碟上的 `<data_dir>/metabolopan/logs/` 下達 7 天，然後在啟動時自動刪除。
若你想擷取先前一次執行的日誌，請在重新開啟應用程式前先到該目錄查找。

<a id="key-references"></a>
## 主要參考文獻（Key references）

- **Brunner–Munzel test.** Brunner, E. & Munzel, U. (2000).
  *The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation.* Biometrical Journal 42 (1): 17–25.
- **Cliff's δ.** Cliff, N. (1993).
  *Dominance statistics: Ordinal analyses to answer ordinal questions.* Psychological Bulletin 114 (3): 494–509.
- **Welch's t-test.** Welch, B. L. (1947).
  *The generalization of "Student's" problem when several different population variances are involved.* Biometrika 34 (1/2): 28–35.
- **BY FDR.** Benjamini, Y. & Yekutieli, D. (2001).
  *The control of the false discovery rate in multiple testing under dependency.* Annals of Statistics 29 (4): 1165–1188.
- **BH FDR.** Benjamini, Y. & Hochberg, Y. (1995).
  *Controlling the false discovery rate: A practical and powerful approach to multiple testing.* Journal of the Royal Statistical Society B 57 (1): 289–300.
- **Hypergeometric ORA.** Khatri, P., Sirota, M. & Butte, A. J. (2012).
  *Ten years of pathway analysis: Current approaches and outstanding challenges.* PLoS Computational Biology 8 (2): e1002375.
- **KEGG.** Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. & Ishiguro-Watanabe, M. (2023).
  *KEGG for taxonomy-based analysis of pathways and genomes.* Nucleic Acids Research 51 (D1): D587–D592.
- **KEGG MODULE.** Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M. & Tanabe, M. (2014).
  *Data, information, knowledge and principle: back to metabolism in KEGG.* Nucleic Acids Research 42 (D1): D199–D205.
- **Quantile normalization.** Bolstad, B. M., Irizarry, R. A., Åstrand, M. & Speed, T. P. (2003).
  *A comparison of normalization methods for high density oligonucleotide array data based on variance and bias.* Bioinformatics 19 (2): 185–193.
- **Probabilistic Quotient Normalization (PQN).** Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. (2006).
  *Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures.
  Application in 1H NMR metabonomics.* Analytical Chemistry 78 (13): 4281–4290.
