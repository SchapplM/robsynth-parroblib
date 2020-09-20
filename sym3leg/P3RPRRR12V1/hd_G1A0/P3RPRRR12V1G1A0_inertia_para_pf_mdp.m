% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR12V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRRR12V1G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:42
% EndTime: 2020-08-06 18:21:43
% DurationCPUTime: 1.42s
% Computational Cost: add. (1722->208), mult. (2193->339), div. (441->14), fcn. (2277->18), ass. (0->169)
t1360 = 0.1e1 / pkin(3);
t1436 = -MDP(9) * t1360 + MDP(4);
t1435 = 2 * MDP(5);
t1434 = 2 * MDP(8);
t1433 = 2 * MDP(12);
t1432 = 2 * MDP(13);
t1359 = pkin(1) + pkin(5);
t1431 = MDP(6) * pkin(1);
t1429 = MDP(10) * t1360;
t1347 = sin(qJ(3,3));
t1320 = t1347 * pkin(3) + qJ(2,3);
t1334 = pkin(6) + t1359;
t1348 = sin(qJ(1,3));
t1354 = cos(qJ(1,3));
t1296 = t1348 * t1320 + t1334 * t1354;
t1299 = -t1320 * t1354 + t1334 * t1348;
t1344 = legFrame(3,3);
t1324 = sin(t1344);
t1327 = cos(t1344);
t1290 = t1296 * t1327 - t1324 * t1299;
t1314 = 0.1e1 / t1320;
t1284 = t1290 * t1314;
t1349 = sin(qJ(3,2));
t1321 = t1349 * pkin(3) + qJ(2,2);
t1350 = sin(qJ(1,2));
t1356 = cos(qJ(1,2));
t1297 = t1350 * t1321 + t1334 * t1356;
t1300 = -t1321 * t1356 + t1334 * t1350;
t1345 = legFrame(2,3);
t1325 = sin(t1345);
t1328 = cos(t1345);
t1291 = t1297 * t1328 - t1325 * t1300;
t1316 = 0.1e1 / t1321;
t1285 = t1291 * t1316;
t1351 = sin(qJ(3,1));
t1322 = -t1351 * pkin(3) - qJ(2,1);
t1352 = sin(qJ(1,1));
t1358 = cos(qJ(1,1));
t1298 = -t1352 * t1322 + t1334 * t1358;
t1301 = t1322 * t1358 + t1334 * t1352;
t1346 = legFrame(1,3);
t1326 = sin(t1346);
t1329 = cos(t1346);
t1292 = t1298 * t1329 - t1326 * t1301;
t1318 = 0.1e1 / t1322;
t1428 = t1292 * t1318;
t1293 = t1324 * t1296 + t1299 * t1327;
t1287 = t1293 * t1314;
t1315 = 0.1e1 / t1320 ^ 2;
t1427 = t1293 * t1315;
t1294 = t1325 * t1297 + t1300 * t1328;
t1288 = t1294 * t1316;
t1317 = 0.1e1 / t1321 ^ 2;
t1426 = t1294 * t1317;
t1295 = t1326 * t1298 + t1301 * t1329;
t1425 = t1295 * t1318;
t1319 = 0.1e1 / t1322 ^ 2;
t1424 = t1295 * t1319;
t1308 = -t1324 * t1348 + t1327 * t1354;
t1302 = t1308 ^ 2;
t1423 = t1302 * t1315;
t1309 = -t1325 * t1350 + t1328 * t1356;
t1303 = t1309 ^ 2;
t1422 = t1303 * t1317;
t1310 = -t1326 * t1352 + t1329 * t1358;
t1304 = t1310 ^ 2;
t1421 = t1304 * t1319;
t1311 = t1324 * t1354 + t1327 * t1348;
t1305 = t1311 ^ 2;
t1420 = t1305 * t1315;
t1312 = t1325 * t1356 + t1328 * t1350;
t1306 = t1312 ^ 2;
t1419 = t1306 * t1317;
t1313 = t1326 * t1358 + t1329 * t1352;
t1307 = t1313 ^ 2;
t1418 = t1307 * t1319;
t1417 = t1308 * t1314;
t1416 = t1309 * t1316;
t1415 = t1310 * t1318;
t1414 = t1311 * t1314;
t1413 = t1311 * t1315;
t1412 = t1312 * t1316;
t1411 = t1312 * t1317;
t1410 = t1313 * t1318;
t1409 = t1313 * t1319;
t1408 = t1314 * t1359;
t1353 = cos(qJ(3,3));
t1341 = t1353 ^ 2;
t1407 = t1315 * t1341;
t1406 = t1316 * t1359;
t1355 = cos(qJ(3,2));
t1342 = t1355 ^ 2;
t1405 = t1317 * t1342;
t1404 = t1318 * t1359;
t1357 = cos(qJ(3,1));
t1343 = t1357 ^ 2;
t1403 = t1319 * t1343;
t1402 = 0.1e1 / t1347 * t1353;
t1401 = 0.1e1 / t1349 * t1355;
t1400 = 0.1e1 / t1351 * t1357;
t1399 = pkin(1) * t1417;
t1398 = pkin(1) * t1416;
t1397 = pkin(1) * t1415;
t1396 = pkin(1) * t1414;
t1395 = pkin(1) * t1412;
t1394 = pkin(1) * t1410;
t1393 = qJ(2,3) * t1423;
t1392 = qJ(2,2) * t1422;
t1391 = qJ(2,1) * t1421;
t1390 = qJ(2,3) * t1420;
t1389 = qJ(2,2) * t1419;
t1388 = qJ(2,1) * t1418;
t1387 = t1308 * t1413;
t1386 = t1309 * t1411;
t1385 = t1310 * t1409;
t1384 = t1314 * t1402;
t1383 = t1315 * t1347 * t1353;
t1382 = t1316 * t1401;
t1381 = t1317 * t1349 * t1355;
t1380 = t1318 * t1400;
t1379 = t1319 * t1351 * t1357;
t1364 = t1308 * t1384 + t1309 * t1382 - t1310 * t1380;
t1369 = -t1415 + t1416 + t1417;
t1378 = t1436 * t1364 + t1369 * t1429;
t1363 = t1311 * t1384 + t1312 * t1382 - t1313 * t1380;
t1368 = -t1410 + t1412 + t1414;
t1377 = t1436 * t1363 + t1368 * t1429;
t1376 = qJ(2,3) * t1387;
t1375 = qJ(2,2) * t1386;
t1374 = qJ(2,1) * t1385;
t1373 = t1353 * t1387;
t1372 = t1355 * t1386;
t1371 = t1357 * t1385;
t1370 = (-t1347 * t1373 - t1349 * t1372 - t1351 * t1371) * t1434 + (t1347 * t1376 + t1349 * t1375 + t1351 * t1374) * t1433 + (qJ(2,1) * t1371 + qJ(2,2) * t1372 + qJ(2,3) * t1373) * t1432 + (t1341 * t1387 + t1342 * t1386 + t1343 * t1385) * MDP(7) + (t1374 + t1375 + t1376) * t1435 + (t1385 + t1386 + t1387) * MDP(1);
t1336 = 0.1e1 / t1347 ^ 2;
t1338 = 0.1e1 / t1349 ^ 2;
t1340 = 0.1e1 / t1351 ^ 2;
t1367 = t1341 * t1336 + t1342 * t1338 + t1343 * t1340;
t1366 = t1290 * t1384 + t1291 * t1382 - t1292 * t1380;
t1365 = t1293 * t1384 + t1294 * t1382 - t1295 * t1380;
t1362 = pkin(1) ^ 2;
t1332 = qJ(2,1) ^ 2 + t1362;
t1331 = qJ(2,2) ^ 2 + t1362;
t1330 = qJ(2,3) ^ 2 + t1362;
t1281 = -t1425 + 0.2e1 * t1394;
t1280 = -t1425 + t1394;
t1279 = t1288 - 0.2e1 * t1395;
t1278 = t1288 - t1395;
t1277 = t1287 - 0.2e1 * t1396;
t1276 = t1287 - t1396;
t1275 = -t1428 + 0.2e1 * t1397;
t1274 = -t1428 + t1397;
t1273 = t1285 - 0.2e1 * t1398;
t1272 = t1285 - t1398;
t1271 = t1284 - 0.2e1 * t1399;
t1270 = t1284 - t1399;
t1269 = -t1313 * t1404 + t1425;
t1268 = t1312 * t1406 - t1288;
t1267 = t1311 * t1408 - t1287;
t1266 = -t1310 * t1404 + t1428;
t1265 = t1309 * t1406 - t1285;
t1264 = t1308 * t1408 - t1284;
t1263 = (-pkin(1) * t1295 + t1313 * t1332) * t1318;
t1262 = (-pkin(1) * t1292 + t1310 * t1332) * t1318;
t1261 = (-pkin(1) * t1294 + t1312 * t1331) * t1316;
t1260 = (-pkin(1) * t1291 + t1309 * t1331) * t1316;
t1259 = (-pkin(1) * t1293 + t1311 * t1330) * t1314;
t1258 = (-pkin(1) * t1290 + t1308 * t1330) * t1314;
t1 = [(t1421 + t1422 + t1423) * MDP(1) + ((-t1275 * t1318 + t1292 * t1319) * t1310 + (t1273 * t1316 + t1291 * t1317) * t1309 + (t1271 * t1314 + t1290 * t1315) * t1308) * MDP(4) + (-(-t1262 * t1310 + t1274 * t1292) * t1318 + (t1260 * t1309 + t1272 * t1291) * t1316 + (t1258 * t1308 + t1270 * t1290) * t1314) * MDP(6) + (t1302 * t1407 + t1303 * t1405 + t1304 * t1403) * MDP(7) + MDP(14) + (t1391 + t1392 + t1393) * t1435 + (-t1302 * t1383 - t1303 * t1381 - t1304 * t1379) * t1434 + (t1347 * t1393 + t1349 * t1392 + t1351 * t1391) * t1433 + (t1353 * t1393 + t1355 * t1392 + t1357 * t1391) * t1432; (t1277 * t1417 + t1279 * t1416 - t1281 * t1415 + t1290 * t1413 + t1291 * t1411 + t1292 * t1409) * MDP(4) + (-(-t1263 * t1310 + t1280 * t1292) * t1318 + (t1261 * t1309 + t1278 * t1291) * t1316 + (t1259 * t1308 + t1276 * t1290) * t1314) * MDP(6) + t1370; t1366 * MDP(6) - t1364 * t1431 + (-t1366 * MDP(12) + (t1285 - t1428 + t1284) * MDP(13) + (MDP(12) * t1364 - MDP(13) * t1369) * t1359) * t1360 + t1378; (t1271 * t1414 + t1273 * t1412 - t1275 * t1410 + t1308 * t1427 + t1309 * t1426 + t1310 * t1424) * MDP(4) + (-(-t1262 * t1313 + t1274 * t1295) * t1318 + (t1260 * t1312 + t1272 * t1294) * t1316 + (t1258 * t1311 + t1270 * t1293) * t1314) * MDP(6) + t1370; (t1418 + t1419 + t1420) * MDP(1) + ((-t1281 * t1318 + t1424) * t1313 + (t1279 * t1316 + t1426) * t1312 + (t1277 * t1314 + t1427) * t1311) * MDP(4) + (-(-t1263 * t1313 + t1280 * t1295) * t1318 + (t1261 * t1312 + t1278 * t1294) * t1316 + (t1259 * t1311 + t1276 * t1293) * t1314) * MDP(6) + (t1305 * t1407 + t1306 * t1405 + t1307 * t1403) * MDP(7) + MDP(14) + (t1388 + t1389 + t1390) * t1435 + (-t1305 * t1383 - t1306 * t1381 - t1307 * t1379) * t1434 + (t1347 * t1390 + t1349 * t1389 + t1351 * t1388) * t1433 + (t1353 * t1390 + t1355 * t1389 + t1357 * t1388) * t1432; t1365 * MDP(6) - t1363 * t1431 + (-t1365 * MDP(12) + (t1288 - t1425 + t1287) * MDP(13) + (MDP(12) * t1363 - MDP(13) * t1368) * t1359) * t1360 + t1377; (t1270 * t1402 + t1272 * t1401 + t1274 * t1400) * MDP(6) + ((t1264 * t1402 + t1265 * t1401 + t1266 * t1400) * MDP(12) + (-t1264 - t1265 - t1266) * MDP(13)) * t1360 + t1378; (t1276 * t1402 + t1278 * t1401 + t1280 * t1400) * MDP(6) + ((t1267 * t1402 + t1268 * t1401 + t1269 * t1400) * MDP(12) + (-t1267 - t1268 - t1269) * MDP(13)) * t1360 + t1377; t1367 * MDP(6) + MDP(14) + (t1336 + t1338 + t1340) * MDP(11) / pkin(3) ^ 2 + 0.2e1 * (-t1367 * MDP(12) + (t1400 + t1401 + t1402) * MDP(13)) * t1360;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
