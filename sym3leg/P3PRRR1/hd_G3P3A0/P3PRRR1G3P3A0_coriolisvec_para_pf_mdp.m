% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G3P3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G3P3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:07:03
% EndTime: 2020-03-09 21:07:05
% DurationCPUTime: 2.24s
% Computational Cost: add. (19020->168), mult. (11613->349), div. (2418->9), fcn. (11529->24), ass. (0->175)
t1411 = pkin(7) + qJ(2,3);
t1402 = qJ(3,3) + t1411;
t1390 = sin(t1402);
t1393 = cos(t1402);
t1423 = xDP(3);
t1414 = legFrame(3,2);
t1405 = sin(t1414);
t1408 = cos(t1414);
t1424 = xDP(2);
t1425 = xDP(1);
t1441 = t1405 * t1424 - t1408 * t1425;
t1354 = t1390 * t1423 + t1441 * t1393;
t1412 = pkin(7) + qJ(2,2);
t1403 = qJ(3,2) + t1412;
t1391 = sin(t1403);
t1394 = cos(t1403);
t1415 = legFrame(2,2);
t1406 = sin(t1415);
t1409 = cos(t1415);
t1440 = t1406 * t1424 - t1409 * t1425;
t1355 = t1391 * t1423 + t1440 * t1394;
t1413 = pkin(7) + qJ(2,1);
t1404 = qJ(3,1) + t1413;
t1392 = sin(t1404);
t1395 = cos(t1404);
t1416 = legFrame(1,2);
t1407 = sin(t1416);
t1410 = cos(t1416);
t1439 = t1407 * t1424 - t1410 * t1425;
t1356 = t1392 * t1423 + t1439 * t1395;
t1396 = sin(t1411);
t1399 = cos(t1411);
t1372 = t1390 * t1399 - t1396 * t1393;
t1533 = 0.1e1 / t1372 ^ 2;
t1542 = t1354 ^ 2 * t1533;
t1397 = sin(t1412);
t1400 = cos(t1412);
t1373 = t1391 * t1400 - t1397 * t1394;
t1532 = 0.1e1 / t1373 ^ 2;
t1541 = t1355 ^ 2 * t1532;
t1398 = sin(t1413);
t1401 = cos(t1413);
t1374 = t1392 * t1401 - t1398 * t1395;
t1531 = 0.1e1 / t1374 ^ 2;
t1540 = t1356 ^ 2 * t1531;
t1536 = 0.1e1 / t1372;
t1539 = t1405 * t1536;
t1535 = 0.1e1 / t1373;
t1538 = t1406 * t1535;
t1534 = 0.1e1 / t1374;
t1537 = t1407 * t1534;
t1386 = pkin(2) * t1401 + pkin(3) * t1395;
t1482 = t1534 * t1386 * t1410;
t1385 = pkin(2) * t1400 + pkin(3) * t1394;
t1484 = t1535 * t1385 * t1409;
t1384 = pkin(2) * t1399 + pkin(3) * t1393;
t1486 = t1536 * t1384 * t1408;
t1428 = 0.1e1 / pkin(2);
t1427 = 0.1e1 / pkin(3);
t1348 = -(-t1396 * t1423 - t1399 * t1441) * pkin(2) + t1354 * pkin(3);
t1527 = t1348 * t1536;
t1497 = t1427 * t1527;
t1524 = t1354 * t1536;
t1343 = (t1497 - t1524) * t1428;
t1426 = pkin(3) ^ 2;
t1444 = t1390 * t1396 + t1393 * t1399;
t1502 = 0.2e1 * pkin(3) * t1428;
t1530 = (-t1343 * t1426 + (t1524 - t1444 * (-t1524 + t1497 / 0.2e1) * t1502) * pkin(2)) * t1427;
t1349 = -(-t1397 * t1423 - t1400 * t1440) * pkin(2) + t1355 * pkin(3);
t1526 = t1349 * t1535;
t1495 = t1427 * t1526;
t1522 = t1355 * t1535;
t1345 = (t1495 - t1522) * t1428;
t1443 = t1391 * t1397 + t1394 * t1400;
t1529 = (-t1345 * t1426 + (t1522 - t1443 * (-t1522 + t1495 / 0.2e1) * t1502) * pkin(2)) * t1427;
t1350 = -(-t1398 * t1423 - t1401 * t1439) * pkin(2) + t1356 * pkin(3);
t1525 = t1350 * t1534;
t1493 = t1427 * t1525;
t1520 = t1356 * t1534;
t1347 = (t1493 - t1520) * t1428;
t1442 = t1392 * t1398 + t1395 * t1401;
t1528 = (-t1347 * t1426 + (t1520 - t1442 * (-t1520 + t1493 / 0.2e1) * t1502) * pkin(2)) * t1427;
t1523 = t1354 * t1533;
t1521 = t1355 * t1532;
t1519 = t1356 * t1531;
t1518 = t1536 * (pkin(2) * t1396 + pkin(3) * t1390);
t1517 = t1535 * (pkin(2) * t1397 + pkin(3) * t1391);
t1516 = t1534 * (pkin(2) * t1398 + pkin(3) * t1392);
t1515 = t1536 * t1390;
t1514 = t1535 * t1391;
t1513 = t1534 * t1392;
t1506 = t1393 * t1408;
t1505 = t1394 * t1409;
t1504 = t1395 * t1410;
t1429 = 0.1e1 / pkin(2) ^ 2;
t1503 = t1427 * t1429;
t1333 = pkin(3) * t1343 - t1444 * t1524;
t1498 = t1536 * t1527;
t1476 = t1343 * t1498;
t1327 = (-t1333 * t1523 + t1476) * t1429;
t1501 = t1327 * t1518;
t1334 = pkin(3) * t1345 - t1443 * t1522;
t1496 = t1535 * t1526;
t1474 = t1345 * t1496;
t1328 = (-t1334 * t1521 + t1474) * t1429;
t1500 = t1328 * t1517;
t1335 = pkin(3) * t1347 - t1442 * t1520;
t1494 = t1534 * t1525;
t1472 = t1347 * t1494;
t1329 = (-t1335 * t1519 + t1472) * t1429;
t1499 = t1329 * t1516;
t1492 = t1384 * t1539;
t1491 = t1536 * t1506;
t1490 = t1385 * t1538;
t1489 = t1535 * t1505;
t1488 = t1386 * t1537;
t1487 = t1534 * t1504;
t1485 = t1393 * t1539;
t1483 = t1394 * t1538;
t1481 = t1395 * t1537;
t1336 = (t1444 * pkin(2) + pkin(3)) * t1476 * t1503;
t1322 = -t1336 + (0.2e1 * t1476 + (-0.2e1 * t1333 - t1530) * t1523) * t1429;
t1480 = t1322 * t1515;
t1337 = (t1443 * pkin(2) + pkin(3)) * t1474 * t1503;
t1324 = -t1337 + (0.2e1 * t1474 + (-0.2e1 * t1334 - t1529) * t1521) * t1429;
t1479 = t1324 * t1514;
t1338 = (t1442 * pkin(2) + pkin(3)) * t1472 * t1503;
t1326 = -t1338 + (0.2e1 * t1472 + (-0.2e1 * t1335 - t1528) * t1519) * t1429;
t1478 = t1326 * t1513;
t1342 = (t1497 - 0.2e1 * t1524) * t1428;
t1477 = t1342 * t1498;
t1344 = (t1495 - 0.2e1 * t1522) * t1428;
t1475 = t1344 * t1496;
t1346 = (t1493 - 0.2e1 * t1520) * t1428;
t1473 = t1346 * t1494;
t1471 = t1518 * t1542;
t1470 = t1486 * t1542;
t1469 = t1517 * t1541;
t1468 = t1484 * t1541;
t1467 = t1516 * t1540;
t1466 = t1482 * t1540;
t1465 = t1322 * t1485;
t1464 = t1324 * t1483;
t1463 = t1326 * t1481;
t1462 = t1327 * t1492;
t1461 = t1328 * t1490;
t1460 = t1329 * t1488;
t1459 = t1327 * t1486;
t1458 = t1322 * t1491;
t1457 = t1328 * t1484;
t1456 = t1324 * t1489;
t1455 = t1329 * t1482;
t1454 = t1326 * t1487;
t1453 = t1390 * t1477;
t1452 = t1342 * t1348 * t1533 * t1506;
t1451 = t1391 * t1475;
t1450 = t1344 * t1349 * t1532 * t1505;
t1449 = t1392 * t1473;
t1448 = t1346 * t1350 * t1531 * t1504;
t1447 = t1492 * t1542;
t1446 = t1490 * t1541;
t1445 = t1488 * t1540;
t1438 = t1393 * t1405 * t1477;
t1437 = t1394 * t1406 * t1475;
t1436 = t1395 * t1407 * t1473;
t1422 = cos(qJ(3,1));
t1421 = cos(qJ(3,2));
t1420 = cos(qJ(3,3));
t1419 = sin(qJ(3,1));
t1418 = sin(qJ(3,2));
t1417 = sin(qJ(3,3));
t1325 = -t1338 + (t1472 + (-t1335 - t1528) * t1519) * t1429;
t1323 = -t1337 + (t1474 + (-t1334 - t1529) * t1521) * t1429;
t1321 = -t1336 + (t1476 + (-t1333 - t1530) * t1523) * t1429;
t1 = [(t1420 * t1458 + t1421 * t1456 + t1422 * t1454) * MDP(6) + (-t1417 * t1458 - t1418 * t1456 - t1419 * t1454) * MDP(7) + ((t1327 * t1491 + t1328 * t1489 + t1329 * t1487) * MDP(2) + (t1321 * t1491 + t1323 * t1489 + t1325 * t1487) * MDP(5)) * t1428 + ((-t1420 * t1459 - t1421 * t1457 - t1422 * t1455) * MDP(6) + (t1417 * t1459 + t1418 * t1457 + t1419 * t1455) * MDP(7) + ((-t1417 * t1470 - t1418 * t1468 - t1419 * t1466) * MDP(6) + (-t1420 * t1470 - t1421 * t1468 - t1422 * t1466) * MDP(7)) * t1429 + ((-t1321 * t1486 - t1323 * t1484 - t1325 * t1482) * MDP(5) + (-t1417 * t1452 - t1418 * t1450 - t1419 * t1448) * MDP(6) + (-t1420 * t1452 - t1421 * t1450 - t1422 * t1448) * MDP(7)) * t1428) * t1427; (-t1420 * t1465 - t1421 * t1464 - t1422 * t1463) * MDP(6) + (t1417 * t1465 + t1418 * t1464 + t1419 * t1463) * MDP(7) + ((-t1327 * t1485 - t1328 * t1483 - t1329 * t1481) * MDP(2) + (-t1321 * t1485 - t1323 * t1483 - t1325 * t1481) * MDP(5)) * t1428 + ((t1420 * t1462 + t1421 * t1461 + t1422 * t1460) * MDP(6) + (-t1417 * t1462 - t1418 * t1461 - t1419 * t1460) * MDP(7) + ((t1417 * t1447 + t1418 * t1446 + t1419 * t1445) * MDP(6) + (t1420 * t1447 + t1421 * t1446 + t1422 * t1445) * MDP(7)) * t1429 + ((t1321 * t1492 + t1323 * t1490 + t1325 * t1488) * MDP(5) + (t1417 * t1438 + t1418 * t1437 + t1419 * t1436) * MDP(6) + (t1420 * t1438 + t1421 * t1437 + t1422 * t1436) * MDP(7)) * t1428) * t1427; (-t1420 * t1480 - t1421 * t1479 - t1422 * t1478) * MDP(6) + (t1417 * t1480 + t1418 * t1479 + t1419 * t1478) * MDP(7) + ((-t1327 * t1515 - t1328 * t1514 - t1329 * t1513) * MDP(2) + (-t1321 * t1515 - t1323 * t1514 - t1325 * t1513) * MDP(5)) * t1428 + ((t1420 * t1501 + t1421 * t1500 + t1422 * t1499) * MDP(6) + (-t1417 * t1501 - t1418 * t1500 - t1419 * t1499) * MDP(7) + ((t1417 * t1471 + t1418 * t1469 + t1419 * t1467) * MDP(6) + (t1420 * t1471 + t1421 * t1469 + t1422 * t1467) * MDP(7)) * t1429 + ((t1321 * t1518 + t1323 * t1517 + t1325 * t1516) * MDP(5) + (t1417 * t1453 + t1418 * t1451 + t1419 * t1449) * MDP(6) + (t1420 * t1453 + t1421 * t1451 + t1422 * t1449) * MDP(7)) * t1428) * t1427;];
taucX  = t1;
