// --- 設定 ---
numPouchPoints = 24; 
saveDir = getDirectory("CSVの保存先フォルダを選んでください");
csvPath = saveDir + "all_samples_coords.csv";

run("Set Measurements...", "display x y redirect=None decimal=3");

if (!File.exists(csvPath)) {
    f = File.open(csvPath);
    print(f, "sample,type,x,y");
    File.close(f);
}

print("\\Clear");
print("--- Analysis Started ---");

ids = getList("image.titles");

for (k = 0; k < ids.length; k++) {
    selectImage(ids[k]);
    imgName = getTitle(); // ここで定義されるので、これより後で print する必要があります
    print("Processing: " + imgName);

    // 1. 向き確認
    setTool("hand");
    waitForUser("Orientation Check", "Zを動かして向きを確認してください。");
    Dialog.create("Flip Options");
    Dialog.addCheckbox("Flip Vertically", false);
    Dialog.addCheckbox("Flip Horizontally", false);
    Dialog.show();
    if (Dialog.getCheckbox()) run("Flip Vertically");
    if (Dialog.getCheckbox()) run("Flip Horizontally");

    // 2. Pouchの縁 (24点) - スタック画像上で指定
    success = false;
    while (!success) {
        setTool("multipoint");
        run("Select None");
        run("Clear Results"); 
        waitForUser("Pouch Edge", "一番見やすいスライスで24点打ってください。");
        run("Measure"); 
        if (nResults == numPouchPoints) {
            success = true;
        } else {
            showMessage("Error", numPouchPoints + "点ちょうどにしてください。");
        }
    }
    px = newArray(nResults); py = newArray(nResults);
    print("[" + imgName + "] type: pouch_edge");
    for (i = 0; i < nResults; i++) {
        px[i] = getResult("X", i);
        py[i] = getResult("Y", i);
        print(px[i] + ", " + py[i]); // Logにバックアップ
    }

    // 3. Z-Projectの範囲指定
    waitForUser("Z-Range", "開始スライスを表示してOK");
    Stack.getPosition(dummyC, startS, dummyF);
    waitForUser("Z-Range", "終了スライスを表示してOK");
    Stack.getPosition(dummyC, endS, dummyF);
    run("Z Project...", "start=" + startS + " stop=" + endS + " projection=[Max Intensity]");
    projName = getTitle();

    // 4. Xrp1の点
    run("Clear Results");
    setTool("multipoint");
    run("Select None");
    waitForUser("Xrp1 Points", "Xrp1をすべて打ってください。");
    run("Measure");
    ax = newArray(nResults); ay = newArray(nResults);
    print("[" + imgName + "] type: Xrp1_pos");
    for (i = 0; i < nResults; i++) {
        ax[i] = getResult("X", i);
        ay[i] = getResult("Y", i);
        print(ax[i] + ", " + ay[i]);
    }

    // 5. cDcp1の点
    //run("Clear Results");
    //setTool("multipoint");
    //run("Select None");
    //waitForUser("cDcp1 Points", "cDcp1をすべて打ってください。");
    //run("Measure");
    //dx = newArray(nResults); dy = newArray(nResults);
    //print("[" + imgName + "] type: cdcp1_pos");
    //for (i = 0; i < nResults; i++) {
    //    dx[i] = getResult("X", i);
    //    dy[i] = getResult("Y", i);
    //    print(dx[i] + ", " + dy[i]);
    //}

    // --- CSVへ追記 ---
    for (i = 0; i < px.length; i++) File.append(imgName + ",pouch_edge," + px[i] + "," + py[i], csvPath);
    for (j = 0; j < ax.length; j++) File.append(imgName + ",xrp1_pos," + ax[j] + "," + ay[j], csvPath);
    //for (m = 0; m < dx.length; m++) File.append(imgName + ",cdcp1_pos," + dx[m] + "," + dy[m], csvPath);
    
    close(); // 投影図を閉じる
    if (isOpen(imgName)) { selectImage(imgName); close(); } // 元画像を閉じる
    print("--- Done: " + imgName + " ---\n");
}
showMessage("Finish!");