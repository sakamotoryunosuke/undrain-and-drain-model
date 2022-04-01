undrain and drain modelは微動の潮汐変調(phase difference, tidal sensitivity)をシミュレートするコードです。
それぞれのfolderにあるnumericalcalculationはdrain/undrain modelを数値的に解いたコードです。
それらのモデルは1自由度バネブロックモデルです。
このモデルはslip lawを用いた状態依存摩擦則によって支配されています。
加えて、shear zoneで生じるdilatancy/compactionを考慮しました。
本モデルは準静的運動方程式を仮定しました。
approximateexpressionはdrain/undrain modelにおけるphase differenceとtidal sensitivityを近似式から導出するコードです。
