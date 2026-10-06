//フォルダ設定
inputDir = "/Volumes/IODATASSD/2024 Igaki lab/Data/2609/260926 Xrp1 upd, acn-RNAi egg 30 wing size/12x/"; 
outputDir = inputDir + "Maxima_Results/";

// 出力フォルダがない場合に作成
if (!File.exists(outputDir)) File.makeDirectory(outputDir);

// Bio-Formatsのマクロ拡張機能を有効化
run("Bio-Formats Macro Extensions");

// 共通設定
size = 236; // ROIサイズ (236x236) = 0.01mm2

setBatchMode(false); // falseをtrueにすると、画面表示を省略して高速化できる
list = getFileList(inputDir);

for (i = 0; i < list.length; i++) {
    // .lifファイルのみを対象にする
    if (endsWith(list[i], ".lif")) {
        path = inputDir + list[i];
        
        // ファイル内の画像シリーズ数を取得
        Ext.setId(path);
        Ext.getSeriesCount(seriesCount);
        
        for (s = 0; s < seriesCount; s++) {
            // Bio-Formatsでシリーズを個別に開く
            run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + (s+1));
            
            // 処理対象の画像情報を取得
            
            fullTitle = getTitle();
            baseName = replace(list[i], ".lif", ""); // 元のファイル名
            print(fullTitle);
            
            // ROIマネージャーをリセット
            roiManager("reset");
            
            // 画像サイズが異なっていた場合
            // run("Scale...", "x=0.5 y=0.5 z=1.0 width=1920 height=1080 depth=3 interpolation=Bicubic average create");
            // selectWindow(fullTitle);
            // close();
            // fullTitle = getTitle();
            // fullTitle = replace(fullTitle, "-1", "");
			// rename(fullTitle);

            // --- 1. 上下反転の判定 (名前の最後が "l" または全体が "l" で始まる場合) ---
            // 元のコードにあった両方の判定ロジックを統合
            parts = split(fullTitle, "-");
            lastPart = trim(parts[lengthOf(parts)-1]);
            
            if (startsWith(fullTitle, "l") || startsWith(lastPart, "l")) {
                showStatus("Flipping " + fullTitle + "...");
                run("Flip Vertically", "stack");
            }

            // --- 2. Anterior(a) か Posterior(p) かを判定して座標をセット ---
            // タイトルの末尾1文字を取得して判定
            if (endsWith(fullTitle, "a")) {
                xPoints = newArray(1300, 1000, 700);
                yPoints = newArray(550, 550, 550);
            } else if (endsWith(fullTitle, "p")) {
                xPoints = newArray(350, 750, 1150);
                yPoints = newArray(580, 580, 580);
            } else {
                print("Skipped: " + fullTitle + " (末尾が a/p ではありません)");
                close();
                continue; // 次のシリーズへ
            }

            // --- 3. 指定された3箇所のROIを処理 ---
            for (j = 0; j < xPoints.length; j++) {
                selectWindow(fullTitle);
                
                // ROIを作成
                makeRectangle(xPoints[j], yPoints[j], size, size);
                
                // 処理用に複製
                run("Duplicate...", "duplicate");
                dupTitle = getTitle();
                
                // 画像処理プロセス
                run("Gaussian Blur...", "sigma=2");
                run("Stack to RGB");
                run("16-bit");
                
                // Find Maximaの実行（カウントを取得したい場合はLogに出力されます）
                run("Find Maxima...", "prominence=1 light output=Count");
                run("Find Maxima...", "prominence=1 light output=[Point Selection]");
                run("Flatten");
                
                // 保存名の設定（「ファイル名_シリーズ名_ROI番号_Maxima.tif」）
                saveName = baseName + "_" + fullTitle + "_ROI" + j + "_Maxima.tif";
                saveAs("Tiff", outputDir + saveName);
                
                // 複製した画像群を閉じる
                close(); // Flatten後の画像
                selectWindow(dupTitle);
                close(); // Duplicateした画像
                close();
            }

            // 最後に元のシリーズ画像を閉じる
            selectWindow(fullTitle);
            close();
        }
    }
}

setBatchMode(false);
print("Done! 全ての処理が完了しました。");