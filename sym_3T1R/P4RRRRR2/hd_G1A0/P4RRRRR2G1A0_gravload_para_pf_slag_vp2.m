% Calculate Gravitation load for parallel robot
% P4RRRRR2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4RRRRR2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:46
% EndTime: 2020-08-07 17:24:46
% DurationCPUTime: 0.96s
% Computational Cost: add. (1383->222), mult. (1861->302), div. (116->18), fcn. (1222->98), ass. (0->171)
t1524 = -2 * pkin(1);
t1458 = legFrame(4,3);
t1430 = sin(t1458);
t1434 = cos(t1458);
t1379 = -t1430 * g(1) + t1434 * g(2);
t1523 = mrSges(3,1) * t1379;
t1459 = legFrame(3,3);
t1431 = sin(t1459);
t1435 = cos(t1459);
t1380 = -t1431 * g(1) + t1435 * g(2);
t1522 = mrSges(3,1) * t1380;
t1460 = legFrame(2,3);
t1432 = sin(t1460);
t1436 = cos(t1460);
t1381 = -t1432 * g(1) + t1436 * g(2);
t1521 = mrSges(3,1) * t1381;
t1461 = legFrame(1,3);
t1433 = sin(t1461);
t1437 = cos(t1461);
t1382 = -t1433 * g(1) + t1437 * g(2);
t1520 = mrSges(3,1) * t1382;
t1519 = mrSges(3,2) * t1379;
t1518 = mrSges(3,2) * t1380;
t1517 = mrSges(3,2) * t1381;
t1516 = mrSges(3,2) * t1382;
t1383 = t1434 * g(1) + t1430 * g(2);
t1515 = mrSges(3,2) * t1383;
t1384 = t1435 * g(1) + t1431 * g(2);
t1514 = mrSges(3,2) * t1384;
t1385 = t1436 * g(1) + t1432 * g(2);
t1513 = mrSges(3,2) * t1385;
t1386 = t1437 * g(1) + t1433 * g(2);
t1512 = mrSges(3,2) * t1386;
t1511 = qJ(2,1) - qJ(3,1);
t1510 = qJ(2,1) + qJ(3,1);
t1509 = qJ(2,2) - qJ(3,2);
t1508 = qJ(2,2) + qJ(3,2);
t1507 = qJ(2,3) - qJ(3,3);
t1506 = qJ(2,3) + qJ(3,3);
t1423 = qJ(1,1) + t1461;
t1422 = qJ(1,2) + t1460;
t1421 = qJ(1,3) + t1459;
t1363 = mrSges(3,1) * t1383;
t1367 = mrSges(2,1) * t1383;
t1415 = (m(2) + m(3)) * pkin(1) + mrSges(1,1);
t1446 = qJ(1,4) + qJ(2,4);
t1417 = sin(t1446);
t1418 = cos(t1446);
t1492 = qJ(2,4) + qJ(3,4);
t1419 = qJ(1,4) + t1492;
t1493 = qJ(2,4) - qJ(3,4);
t1420 = qJ(1,4) + t1493;
t1464 = mrSges(3,3) - mrSges(2,2);
t1497 = t1464 * t1418;
t1347 = (-t1515 - t1523) * cos(t1420) / 0.2e1 + (t1363 - t1519) * sin(t1420) / 0.2e1 + (t1515 - t1523) * cos(t1419) / 0.2e1 + (t1363 + t1519) * sin(t1419) / 0.2e1 - mrSges(2,1) * t1379 * t1418 + (t1379 * mrSges(1,2) + t1415 * t1383) * sin(qJ(1,4)) - t1383 * t1497 + (-t1379 * t1464 + t1367) * t1417 + (mrSges(1,2) * t1383 - t1415 * t1379) * cos(qJ(1,4));
t1447 = 0.1e1 / sin(qJ(2,4));
t1505 = t1347 * t1447;
t1364 = mrSges(3,1) * t1384;
t1368 = mrSges(2,1) * t1384;
t1455 = qJ(1,3) + qJ(2,3);
t1424 = sin(t1455);
t1427 = cos(t1455);
t1438 = qJ(1,3) + t1506;
t1439 = qJ(1,3) + t1507;
t1496 = t1464 * t1427;
t1348 = (-t1514 - t1522) * cos(t1439) / 0.2e1 + (t1364 - t1518) * sin(t1439) / 0.2e1 + (t1514 - t1522) * cos(t1438) / 0.2e1 + (t1364 + t1518) * sin(t1438) / 0.2e1 - mrSges(2,1) * t1380 * t1427 + (t1380 * mrSges(1,2) + t1415 * t1384) * sin(qJ(1,3)) - t1384 * t1496 + (-t1380 * t1464 + t1368) * t1424 + (mrSges(1,2) * t1384 - t1415 * t1380) * cos(qJ(1,3));
t1449 = 0.1e1 / sin(qJ(2,3));
t1504 = t1348 * t1449;
t1365 = mrSges(3,1) * t1385;
t1369 = mrSges(2,1) * t1385;
t1456 = qJ(2,2) + qJ(1,2);
t1425 = sin(t1456);
t1428 = cos(t1456);
t1440 = qJ(1,2) + t1508;
t1441 = qJ(1,2) + t1509;
t1495 = t1464 * t1428;
t1349 = (-t1513 - t1521) * cos(t1441) / 0.2e1 + (t1365 - t1517) * sin(t1441) / 0.2e1 + (t1513 - t1521) * cos(t1440) / 0.2e1 + (t1365 + t1517) * sin(t1440) / 0.2e1 - mrSges(2,1) * t1381 * t1428 + (t1381 * mrSges(1,2) + t1415 * t1385) * sin(qJ(1,2)) - t1385 * t1495 + (-t1381 * t1464 + t1369) * t1425 + (mrSges(1,2) * t1385 - t1415 * t1381) * cos(qJ(1,2));
t1450 = 0.1e1 / sin(qJ(2,2));
t1503 = t1349 * t1450;
t1366 = mrSges(3,1) * t1386;
t1370 = mrSges(2,1) * t1386;
t1457 = qJ(1,1) + qJ(2,1);
t1426 = sin(t1457);
t1429 = cos(t1457);
t1442 = qJ(1,1) + t1510;
t1443 = qJ(1,1) + t1511;
t1494 = t1464 * t1429;
t1350 = (-t1512 - t1520) * cos(t1443) / 0.2e1 + (t1366 - t1516) * sin(t1443) / 0.2e1 + (t1512 - t1520) * cos(t1442) / 0.2e1 + (t1366 + t1516) * sin(t1442) / 0.2e1 - mrSges(2,1) * t1382 * t1429 + (t1382 * mrSges(1,2) + t1415 * t1386) * sin(qJ(1,1)) - t1386 * t1494 + (-t1382 * t1464 + t1370) * t1426 + (mrSges(1,2) * t1386 - t1415 * t1382) * cos(qJ(1,1));
t1451 = 0.1e1 / sin(qJ(2,1));
t1502 = t1350 * t1451;
t1462 = sin(qJ(3,4));
t1463 = cos(qJ(3,4));
t1487 = t1463 * mrSges(3,1) - t1462 * mrSges(3,2);
t1351 = t1367 * t1417 + (t1487 * t1417 - t1497) * t1383 + ((-mrSges(2,1) - t1487) * t1418 - t1464 * t1417) * t1379;
t1501 = t1351 / (sin(t1492) + sin(t1493));
t1465 = sin(qJ(3,3));
t1468 = cos(qJ(3,3));
t1486 = t1468 * mrSges(3,1) - t1465 * mrSges(3,2);
t1352 = t1368 * t1424 + (t1486 * t1424 - t1496) * t1384 + ((-mrSges(2,1) - t1486) * t1427 - t1464 * t1424) * t1380;
t1500 = t1352 / (sin(t1506) + sin(t1507));
t1466 = sin(qJ(3,2));
t1469 = cos(qJ(3,2));
t1485 = t1469 * mrSges(3,1) - t1466 * mrSges(3,2);
t1353 = t1369 * t1425 + (t1485 * t1425 - t1495) * t1385 + ((-mrSges(2,1) - t1485) * t1428 - t1464 * t1425) * t1381;
t1499 = t1353 / (sin(t1508) + sin(t1509));
t1467 = sin(qJ(3,1));
t1470 = cos(qJ(3,1));
t1484 = t1470 * mrSges(3,1) - t1467 * mrSges(3,2);
t1354 = t1370 * t1426 + (t1484 * t1426 - t1494) * t1386 + ((-mrSges(2,1) - t1484) * t1429 - t1464 * t1426) * t1382;
t1498 = t1354 / (sin(t1510) + sin(t1511));
t1416 = qJ(1,4) + t1458;
t1414 = qJ(2,1) + t1423;
t1413 = qJ(2,2) + t1422;
t1412 = qJ(2,3) + t1421;
t1483 = 1 / pkin(1);
t1491 = t1447 * t1462 * t1483;
t1490 = t1449 * t1465 * t1483;
t1489 = t1450 * t1466 * t1483;
t1488 = t1451 * t1467 * t1483;
t1411 = qJ(2,4) + t1416;
t1482 = 0.1e1 / pkin(2);
t1481 = koppelP(1,1);
t1480 = koppelP(2,1);
t1479 = koppelP(3,1);
t1478 = koppelP(4,1);
t1477 = koppelP(1,2);
t1476 = koppelP(2,2);
t1475 = koppelP(3,2);
t1474 = koppelP(4,2);
t1473 = mrSges(4,1);
t1472 = mrSges(4,2);
t1471 = xP(4);
t1454 = 0.1e1 / t1470;
t1453 = 0.1e1 / t1469;
t1452 = 0.1e1 / t1468;
t1448 = 0.1e1 / t1463;
t1445 = cos(t1471);
t1444 = sin(t1471);
t1410 = -qJ(3,1) + t1414;
t1409 = qJ(3,1) + t1414;
t1408 = -qJ(3,2) + t1413;
t1407 = qJ(3,2) + t1413;
t1406 = -qJ(3,3) + t1412;
t1405 = qJ(3,3) + t1412;
t1404 = cos(t1414);
t1403 = cos(t1413);
t1402 = cos(t1412);
t1401 = sin(t1414);
t1400 = sin(t1413);
t1399 = sin(t1412);
t1398 = -qJ(3,4) + t1411;
t1397 = qJ(3,4) + t1411;
t1396 = cos(t1411);
t1395 = sin(t1411);
t1378 = -t1444 * t1477 + t1445 * t1481;
t1377 = -t1444 * t1476 + t1445 * t1480;
t1376 = -t1444 * t1475 + t1445 * t1479;
t1375 = -t1444 * t1474 + t1445 * t1478;
t1374 = -t1444 * t1481 - t1445 * t1477;
t1373 = -t1444 * t1480 - t1445 * t1476;
t1372 = -t1444 * t1479 - t1445 * t1475;
t1371 = -t1444 * t1478 - t1445 * t1474;
t1362 = cos(t1423) * t1524 + (-cos(t1409) - cos(t1410)) * pkin(2);
t1361 = cos(t1422) * t1524 + (-cos(t1407) - cos(t1408)) * pkin(2);
t1360 = cos(t1421) * t1524 + (-cos(t1405) - cos(t1406)) * pkin(2);
t1359 = sin(t1423) * t1524 + (-sin(t1410) - sin(t1409)) * pkin(2);
t1358 = sin(t1422) * t1524 + (-sin(t1408) - sin(t1407)) * pkin(2);
t1357 = sin(t1421) * t1524 + (-sin(t1406) - sin(t1405)) * pkin(2);
t1356 = cos(t1416) * t1524 + (-cos(t1397) - cos(t1398)) * pkin(2);
t1355 = sin(t1416) * t1524 + (-sin(t1398) - sin(t1397)) * pkin(2);
t1 = [-g(1) * m(4) + (t1396 * t1505 + t1402 * t1504 + t1403 * t1503 + t1404 * t1502 + (t1356 * t1501 + t1360 * t1500 + t1361 * t1499 + t1362 * t1498) * t1482) * t1483; -g(2) * m(4) + (t1395 * t1505 + t1399 * t1504 + t1400 * t1503 + t1401 * t1502 + (t1355 * t1501 + t1357 * t1500 + t1358 * t1499 + t1359 * t1498) * t1482) * t1483; t1448 * t1347 * t1491 + t1452 * t1348 * t1490 + t1453 * t1349 * t1489 + t1454 * t1350 * t1488 - g(3) * m(4) + (-(cos(qJ(2,1)) * pkin(1) + t1470 * pkin(2)) / t1470 ^ 2 * t1354 * t1488 + t1454 * (-g(3) * t1484 + (t1382 * t1426 + t1386 * t1429) * (mrSges(3,1) * t1467 + mrSges(3,2) * t1470)) - (cos(qJ(2,2)) * pkin(1) + t1469 * pkin(2)) / t1469 ^ 2 * t1353 * t1489 + t1453 * (-g(3) * t1485 + (t1381 * t1425 + t1385 * t1428) * (mrSges(3,1) * t1466 + mrSges(3,2) * t1469)) - (cos(qJ(2,3)) * pkin(1) + t1468 * pkin(2)) / t1468 ^ 2 * t1352 * t1490 + t1452 * (-g(3) * t1486 + (t1380 * t1424 + t1384 * t1427) * (mrSges(3,1) * t1465 + mrSges(3,2) * t1468)) - (cos(qJ(2,4)) * pkin(1) + t1463 * pkin(2)) / t1463 ^ 2 * t1351 * t1491 + t1448 * (-g(3) * t1487 + (t1379 * t1417 + t1383 * t1418) * (mrSges(3,1) * t1462 + mrSges(3,2) * t1463))) * t1482; -(-g(1) * t1473 - g(2) * t1472) * t1444 + t1445 * (g(1) * t1472 - g(2) * t1473) + ((t1374 * t1404 + t1378 * t1401) * t1502 + (t1373 * t1403 + t1377 * t1400) * t1503 + (t1372 * t1402 + t1376 * t1399) * t1504 + (t1371 * t1396 + t1375 * t1395) * t1505 + ((t1359 * t1378 + t1362 * t1374) * t1498 + (t1358 * t1377 + t1361 * t1373) * t1499 + (t1357 * t1376 + t1360 * t1372) * t1500 + (t1355 * t1375 + t1356 * t1371) * t1501) * t1482) * t1483;];
taugX  = t1;
