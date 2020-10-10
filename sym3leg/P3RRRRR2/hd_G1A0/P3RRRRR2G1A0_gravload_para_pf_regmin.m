% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:43
% EndTime: 2020-03-09 21:04:44
% DurationCPUTime: 0.65s
% Computational Cost: add. (903->144), mult. (855->240), div. (255->14), fcn. (927->60), ass. (0->133)
t1394 = legFrame(3,3);
t1376 = sin(t1394);
t1455 = cos(t1394);
t1337 = t1376 * g(1) - t1455 * g(2);
t1340 = t1455 * g(1) + t1376 * g(2);
t1391 = qJ(1,3) + qJ(2,3);
t1370 = sin(t1391);
t1373 = cos(t1391);
t1319 = t1337 * t1373 + t1340 * t1370;
t1397 = sin(qJ(3,3));
t1459 = t1319 * t1397 ^ 2;
t1395 = legFrame(2,3);
t1377 = sin(t1395);
t1454 = cos(t1395);
t1338 = t1377 * g(1) - t1454 * g(2);
t1341 = t1454 * g(1) + t1377 * g(2);
t1392 = qJ(1,2) + qJ(2,2);
t1371 = sin(t1392);
t1374 = cos(t1392);
t1320 = t1338 * t1374 + t1341 * t1371;
t1399 = sin(qJ(3,2));
t1458 = t1320 * t1399 ^ 2;
t1396 = legFrame(1,3);
t1378 = sin(t1396);
t1453 = cos(t1396);
t1339 = t1378 * g(1) - t1453 * g(2);
t1342 = t1453 * g(1) + t1378 * g(2);
t1393 = qJ(1,1) + qJ(2,1);
t1372 = sin(t1393);
t1375 = cos(t1393);
t1321 = t1339 * t1375 + t1342 * t1372;
t1401 = sin(qJ(3,1));
t1457 = t1321 * t1401 ^ 2;
t1456 = -2 * pkin(1);
t1369 = qJ(1,1) + t1396;
t1368 = qJ(1,2) + t1395;
t1367 = qJ(1,3) + t1394;
t1403 = cos(qJ(3,3));
t1452 = t1319 * t1403;
t1405 = cos(qJ(3,2));
t1451 = t1320 * t1405;
t1407 = cos(qJ(3,1));
t1450 = t1321 * t1407;
t1361 = qJ(2,3) + t1367;
t1355 = qJ(3,3) + t1361;
t1356 = -qJ(3,3) + t1361;
t1331 = sin(t1367) * t1456 + (-sin(t1356) - sin(t1355)) * pkin(2);
t1343 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t1449 = t1331 * t1343;
t1362 = qJ(2,2) + t1368;
t1357 = qJ(3,2) + t1362;
t1358 = -qJ(3,2) + t1362;
t1332 = sin(t1368) * t1456 + (-sin(t1358) - sin(t1357)) * pkin(2);
t1344 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t1448 = t1332 * t1344;
t1363 = qJ(2,1) + t1369;
t1359 = qJ(3,1) + t1363;
t1360 = -qJ(3,1) + t1363;
t1333 = sin(t1369) * t1456 + (-sin(t1360) - sin(t1359)) * pkin(2);
t1345 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t1447 = t1333 * t1345;
t1334 = cos(t1367) * t1456 + (-cos(t1356) - cos(t1355)) * pkin(2);
t1446 = t1334 * t1343;
t1335 = cos(t1368) * t1456 + (-cos(t1358) - cos(t1357)) * pkin(2);
t1445 = t1335 * t1344;
t1336 = cos(t1369) * t1456 + (-cos(t1360) - cos(t1359)) * pkin(2);
t1444 = t1336 * t1345;
t1349 = sin(t1361);
t1380 = 0.1e1 / sin(qJ(2,3));
t1443 = t1349 * t1380;
t1350 = sin(t1362);
t1382 = 0.1e1 / sin(qJ(2,2));
t1442 = t1350 * t1382;
t1351 = sin(t1363);
t1384 = 0.1e1 / sin(qJ(2,1));
t1441 = t1351 * t1384;
t1352 = cos(t1361);
t1440 = t1352 * t1380;
t1353 = cos(t1362);
t1439 = t1353 * t1382;
t1354 = cos(t1363);
t1438 = t1354 * t1384;
t1437 = t1380 * t1397;
t1436 = t1382 * t1399;
t1435 = t1384 * t1401;
t1434 = t1319 * t1343 * t1397;
t1433 = t1343 * t1452;
t1432 = t1380 * t1452;
t1431 = t1320 * t1344 * t1399;
t1430 = t1344 * t1451;
t1429 = t1382 * t1451;
t1428 = t1321 * t1345 * t1401;
t1427 = t1345 * t1450;
t1426 = t1384 * t1450;
t1346 = cos(qJ(2,3)) * pkin(1) + t1403 * pkin(2);
t1425 = t1346 * t1380 / t1403 ^ 2;
t1347 = cos(qJ(2,2)) * pkin(1) + t1405 * pkin(2);
t1424 = t1347 * t1382 / t1405 ^ 2;
t1348 = cos(qJ(2,1)) * pkin(1) + t1407 * pkin(2);
t1423 = t1348 * t1384 / t1407 ^ 2;
t1385 = 0.1e1 / t1403;
t1422 = t1385 * t1437;
t1387 = 0.1e1 / t1405;
t1421 = t1387 * t1436;
t1389 = 0.1e1 / t1407;
t1420 = t1389 * t1435;
t1419 = t1319 * t1437;
t1418 = t1320 * t1436;
t1417 = t1321 * t1435;
t1416 = t1385 * t1419;
t1415 = t1387 * t1418;
t1414 = t1389 * t1417;
t1413 = t1397 * t1425;
t1412 = t1399 * t1424;
t1411 = t1401 * t1423;
t1322 = -t1337 * t1370 + t1340 * t1373;
t1323 = -t1338 * t1371 + t1341 * t1374;
t1324 = -t1339 * t1372 + t1342 * t1375;
t1410 = 1 / pkin(1);
t1409 = 0.1e1 / pkin(2);
t1408 = cos(qJ(1,1));
t1406 = cos(qJ(1,2));
t1404 = cos(qJ(1,3));
t1402 = sin(qJ(1,1));
t1400 = sin(qJ(1,2));
t1398 = sin(qJ(1,3));
t1330 = -t1339 * t1402 + t1342 * t1408;
t1329 = -t1338 * t1400 + t1341 * t1406;
t1328 = -t1337 * t1398 + t1340 * t1404;
t1327 = t1339 * t1408 + t1342 * t1402;
t1326 = t1338 * t1406 + t1341 * t1400;
t1325 = t1337 * t1404 + t1340 * t1398;
t1 = [0, (t1325 * t1440 + t1326 * t1439 + t1327 * t1438) * t1410, (t1328 * t1440 + t1329 * t1439 + t1330 * t1438) * t1410, 0, (t1319 * t1440 + t1320 * t1439 + t1321 * t1438 + (t1319 * t1446 + t1320 * t1445 + t1321 * t1444) * t1409) * t1410, (t1322 * t1440 + t1323 * t1439 + t1324 * t1438 + (t1322 * t1446 + t1323 * t1445 + t1324 * t1444) * t1409) * t1410, 0, 0, 0, 0, 0, (t1352 * t1432 + t1353 * t1429 + t1354 * t1426 + (t1334 * t1433 + t1335 * t1430 + t1336 * t1427) * t1409) * t1410, (-t1352 * t1419 - t1353 * t1418 - t1354 * t1417 + (-t1334 * t1434 - t1335 * t1431 - t1336 * t1428) * t1409) * t1410, -g(1); 0, (t1325 * t1443 + t1326 * t1442 + t1327 * t1441) * t1410, (t1328 * t1443 + t1329 * t1442 + t1330 * t1441) * t1410, 0, (t1319 * t1443 + t1320 * t1442 + t1321 * t1441 + (t1319 * t1449 + t1320 * t1448 + t1321 * t1447) * t1409) * t1410, (t1322 * t1443 + t1323 * t1442 + t1324 * t1441 + (t1322 * t1449 + t1323 * t1448 + t1324 * t1447) * t1409) * t1410, 0, 0, 0, 0, 0, (t1349 * t1432 + t1350 * t1429 + t1351 * t1426 + (t1331 * t1433 + t1332 * t1430 + t1333 * t1427) * t1409) * t1410, (-t1349 * t1419 - t1350 * t1418 - t1351 * t1417 + (-t1331 * t1434 - t1332 * t1431 - t1333 * t1428) * t1409) * t1410, -g(2); 0, (t1325 * t1422 + t1326 * t1421 + t1327 * t1420) * t1410, (t1328 * t1422 + t1329 * t1421 + t1330 * t1420) * t1410, 0, (t1416 + t1415 + t1414 + (-t1319 * t1413 - t1320 * t1412 - t1321 * t1411) * t1409) * t1410, (t1322 * t1422 + t1323 * t1421 + t1324 * t1420 + (-t1322 * t1413 - t1323 * t1412 - t1324 * t1411) * t1409) * t1410, 0, 0, 0, 0, 0, (t1417 + t1418 + t1419) * t1410 + (t1389 * (-g(3) * t1407 + t1324 * t1401) + t1387 * (-g(3) * t1405 + t1323 * t1399) + t1385 * (-g(3) * t1403 + t1322 * t1397) + (-t1346 * t1416 - t1347 * t1415 - t1348 * t1414) * t1410) * t1409, (-t1385 * t1380 * t1459 - t1387 * t1382 * t1458 - t1389 * t1384 * t1457) * t1410 + (t1389 * (g(3) * t1401 + t1324 * t1407) + t1387 * (g(3) * t1399 + t1323 * t1405) + t1385 * (g(3) * t1397 + t1322 * t1403) + (t1423 * t1457 + t1424 * t1458 + t1425 * t1459) * t1410) * t1409, -g(3);];
tau_reg  = t1;
