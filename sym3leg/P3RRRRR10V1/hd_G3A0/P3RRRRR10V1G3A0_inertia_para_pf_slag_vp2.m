% Calculate inertia matrix for parallel robot
% P3RRRRR10V1G3A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 23:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 23:01:48
% EndTime: 2020-08-06 23:01:57
% DurationCPUTime: 9.08s
% Computational Cost: add. (11325->596), mult. (27783->1072), div. (648->10), fcn. (21699->26), ass. (0->431)
t1416 = mrSges(3,3) - mrSges(2,2);
t1717 = t1416 * pkin(5) + Ifges(2,6);
t1420 = sin(qJ(3,3));
t1429 = cos(qJ(3,3));
t1331 = Ifges(3,5) * t1420 + Ifges(3,6) * t1429;
t1404 = t1429 ^ 2;
t1415 = Ifges(3,2) - Ifges(3,1);
t1421 = sin(qJ(2,3));
t1430 = cos(qJ(2,3));
t1608 = t1420 * t1429;
t1687 = Ifges(2,5) - Ifges(3,4);
t1691 = 2 * Ifges(3,4);
t1337 = t1420 * mrSges(3,2) - mrSges(2,1);
t1707 = -mrSges(3,1) * t1429 + t1337;
t1716 = t1430 * (-t1331 + t1717) + (pkin(5) * t1707 + t1404 * t1691 - t1415 * t1608 + t1687) * t1421;
t1423 = sin(qJ(3,2));
t1432 = cos(qJ(3,2));
t1332 = Ifges(3,5) * t1423 + Ifges(3,6) * t1432;
t1407 = t1432 ^ 2;
t1424 = sin(qJ(2,2));
t1433 = cos(qJ(2,2));
t1604 = t1423 * t1432;
t1338 = t1423 * mrSges(3,2) - mrSges(2,1);
t1706 = -mrSges(3,1) * t1432 + t1338;
t1715 = t1433 * (-t1332 + t1717) + (pkin(5) * t1706 + t1407 * t1691 - t1415 * t1604 + t1687) * t1424;
t1426 = sin(qJ(3,1));
t1435 = cos(qJ(3,1));
t1333 = Ifges(3,5) * t1426 + Ifges(3,6) * t1435;
t1410 = t1435 ^ 2;
t1427 = sin(qJ(2,1));
t1436 = cos(qJ(2,1));
t1600 = t1426 * t1435;
t1339 = t1426 * mrSges(3,2) - mrSges(2,1);
t1705 = -mrSges(3,1) * t1435 + t1339;
t1714 = t1436 * (-t1333 + t1717) + (pkin(5) * t1705 + t1410 * t1691 - t1415 * t1600 + t1687) * t1427;
t1414 = cos(pkin(3));
t1400 = t1414 ^ 2;
t1686 = -t1400 + 0.1e1;
t1385 = t1421 * pkin(6);
t1349 = t1385 + pkin(1);
t1386 = t1424 * pkin(6);
t1353 = t1386 + pkin(1);
t1387 = t1427 * pkin(6);
t1357 = t1387 + pkin(1);
t1422 = sin(qJ(1,3));
t1667 = t1420 * pkin(5);
t1346 = pkin(2) + t1667;
t1683 = pkin(2) * t1404;
t1513 = -t1346 + 0.2e1 * t1683;
t1713 = t1422 * t1513;
t1425 = sin(qJ(1,2));
t1665 = t1423 * pkin(5);
t1350 = pkin(2) + t1665;
t1682 = pkin(2) * t1407;
t1512 = -t1350 + 0.2e1 * t1682;
t1712 = t1425 * t1512;
t1428 = sin(qJ(1,1));
t1663 = t1426 * pkin(5);
t1354 = pkin(2) + t1663;
t1681 = pkin(2) * t1410;
t1511 = -t1354 + 0.2e1 * t1681;
t1711 = t1428 * t1511;
t1555 = Ifges(3,4) * t1608;
t1334 = 0.2e1 * t1555;
t1612 = t1415 * t1404;
t1688 = Ifges(3,1) + Ifges(2,3);
t1710 = t1334 + t1612 + t1688;
t1554 = Ifges(3,4) * t1604;
t1335 = 0.2e1 * t1554;
t1611 = t1415 * t1407;
t1709 = t1335 + t1611 + t1688;
t1553 = Ifges(3,4) * t1600;
t1336 = 0.2e1 * t1553;
t1610 = t1415 * t1410;
t1708 = t1336 + t1610 + t1688;
t1689 = Ifges(2,1) + Ifges(3,2);
t1514 = Ifges(2,2) + Ifges(3,3) - t1689;
t1704 = -t1514 - t1610;
t1703 = -t1514 - t1611;
t1702 = -t1514 - t1612;
t1593 = t1435 * t1436;
t1701 = -pkin(2) * t1593 - t1387;
t1595 = t1432 * t1433;
t1700 = -pkin(2) * t1595 - t1386;
t1597 = t1429 * t1430;
t1699 = -pkin(2) * t1597 - t1385;
t1441 = m(2) + m(3);
t1698 = (-t1441 * pkin(5) - (2 * mrSges(2,3))) * pkin(5);
t1431 = cos(qJ(1,3));
t1413 = sin(pkin(3));
t1670 = pkin(6) * t1413;
t1501 = (t1430 + 0.1e1) * (t1430 - 0.1e1) * t1670;
t1679 = pkin(2) * t1421;
t1564 = t1422 * t1679;
t1619 = t1413 * t1431;
t1694 = (-(pkin(1) * t1619 + t1422 * pkin(5)) * t1421 + t1431 * t1501) * t1420 - t1564;
t1434 = cos(qJ(1,2));
t1500 = (t1433 + 0.1e1) * (t1433 - 0.1e1) * t1670;
t1678 = pkin(2) * t1424;
t1562 = t1425 * t1678;
t1618 = t1413 * t1434;
t1693 = (-(pkin(1) * t1618 + t1425 * pkin(5)) * t1424 + t1434 * t1500) * t1423 - t1562;
t1437 = cos(qJ(1,1));
t1499 = (t1436 + 0.1e1) * (t1436 - 0.1e1) * t1670;
t1677 = pkin(2) * t1427;
t1560 = t1428 * t1677;
t1617 = t1413 * t1437;
t1692 = (-(pkin(1) * t1617 + t1428 * pkin(5)) * t1427 + t1437 * t1499) * t1426 - t1560;
t1685 = pkin(1) * t1413;
t1684 = pkin(1) * t1414;
t1388 = pkin(1) * t1421;
t1389 = pkin(1) * t1424;
t1390 = pkin(1) * t1427;
t1680 = pkin(2) * t1413;
t1676 = pkin(2) * t1429;
t1675 = pkin(2) * t1432;
t1674 = pkin(2) * t1435;
t1673 = pkin(5) * t1430;
t1672 = pkin(5) * t1433;
t1671 = pkin(5) * t1436;
t1668 = t1420 * pkin(2);
t1666 = t1423 * pkin(2);
t1664 = t1426 * pkin(2);
t1662 = t1430 * pkin(6);
t1661 = t1433 * pkin(6);
t1660 = t1436 * pkin(6);
t1656 = Ifges(3,4) * t1420;
t1655 = Ifges(3,4) * t1423;
t1654 = Ifges(3,4) * t1426;
t1653 = Ifges(3,6) * t1420 + Ifges(2,4);
t1652 = Ifges(3,6) * t1423 + Ifges(2,4);
t1651 = Ifges(3,6) * t1426 + Ifges(2,4);
t1406 = t1430 ^ 2;
t1325 = (t1406 - 0.2e1) * t1668 - pkin(5);
t1295 = t1325 * t1422 * t1413 - t1349 * t1431;
t1322 = -pkin(5) + (t1406 - 0.1e1) * t1668;
t1417 = legFrame(3,2);
t1379 = sin(t1417);
t1382 = cos(t1417);
t1606 = t1421 * t1430;
t1519 = t1420 * t1606;
t1468 = t1519 * t1680;
t1451 = t1431 * t1468;
t1569 = t1422 * t1662;
t1498 = t1379 * t1569;
t1504 = t1404 * t1564;
t1607 = t1420 * t1430;
t1526 = t1382 * t1607;
t1572 = t1413 * t1662;
t1596 = t1430 * t1431;
t1622 = t1413 * t1421;
t1244 = ((-t1325 * t1382 - t1498) * t1429 + (-pkin(6) * t1526 + t1379 * t1713) * t1421) * t1400 + (-(t1379 * t1596 - 0.2e1 * t1382 * t1622) * t1683 + (t1295 * t1379 - t1382 * t1572) * t1429 + (-t1382 * t1346 + t1420 * t1498) * t1622) * t1414 - t1379 * t1504 + (t1382 * t1322 + t1379 * t1451) * t1429 - t1694 * t1379 + t1349 * t1526;
t1405 = 0.1e1 / t1429;
t1650 = t1244 * t1405;
t1409 = t1433 ^ 2;
t1326 = (t1409 - 0.2e1) * t1666 - pkin(5);
t1296 = t1326 * t1425 * t1413 - t1353 * t1434;
t1323 = -pkin(5) + (t1409 - 0.1e1) * t1666;
t1418 = legFrame(2,2);
t1380 = sin(t1418);
t1383 = cos(t1418);
t1602 = t1424 * t1433;
t1518 = t1423 * t1602;
t1467 = t1518 * t1680;
t1450 = t1434 * t1467;
t1567 = t1425 * t1661;
t1496 = t1380 * t1567;
t1503 = t1407 * t1562;
t1603 = t1423 * t1433;
t1516 = t1353 * t1603;
t1568 = pkin(6) * t1603;
t1571 = t1413 * t1661;
t1594 = t1433 * t1434;
t1621 = t1413 * t1424;
t1245 = ((-t1326 * t1383 - t1496) * t1432 + (t1380 * t1712 - t1383 * t1568) * t1424) * t1400 + (-(t1380 * t1594 - 0.2e1 * t1383 * t1621) * t1682 + (t1296 * t1380 - t1383 * t1571) * t1432 + (-t1383 * t1350 + t1423 * t1496) * t1621) * t1414 - t1380 * t1503 + (t1383 * t1323 + t1380 * t1450) * t1432 - t1693 * t1380 + t1383 * t1516;
t1408 = 0.1e1 / t1432;
t1649 = t1245 * t1408;
t1412 = t1436 ^ 2;
t1327 = (t1412 - 0.2e1) * t1664 - pkin(5);
t1297 = t1327 * t1428 * t1413 - t1357 * t1437;
t1324 = -pkin(5) + (t1412 - 0.1e1) * t1664;
t1419 = legFrame(1,2);
t1381 = sin(t1419);
t1384 = cos(t1419);
t1598 = t1427 * t1436;
t1517 = t1426 * t1598;
t1466 = t1517 * t1680;
t1449 = t1437 * t1466;
t1565 = t1428 * t1660;
t1494 = t1381 * t1565;
t1502 = t1410 * t1560;
t1599 = t1426 * t1436;
t1515 = t1357 * t1599;
t1566 = pkin(6) * t1599;
t1570 = t1413 * t1660;
t1592 = t1436 * t1437;
t1620 = t1413 * t1427;
t1246 = ((-t1327 * t1384 - t1494) * t1435 + (t1381 * t1711 - t1384 * t1566) * t1427) * t1400 + (-(t1381 * t1592 - 0.2e1 * t1384 * t1620) * t1681 + (t1297 * t1381 - t1384 * t1570) * t1435 + (-t1384 * t1354 + t1426 * t1494) * t1620) * t1414 - t1381 * t1502 + (t1384 * t1324 + t1381 * t1449) * t1435 - t1692 * t1381 + t1384 * t1515;
t1411 = 0.1e1 / t1435;
t1648 = t1246 * t1411;
t1497 = t1382 * t1569;
t1527 = t1379 * t1607;
t1247 = ((-t1325 * t1379 + t1497) * t1429 + (-pkin(6) * t1527 - t1382 * t1713) * t1421) * t1400 + ((0.2e1 * t1379 * t1622 + t1382 * t1596) * t1683 + (-t1295 * t1382 - t1379 * t1572) * t1429 - (t1379 * t1346 + t1420 * t1497) * t1622) * t1414 + t1382 * t1504 + (t1379 * t1322 - t1382 * t1451) * t1429 + t1694 * t1382 + t1349 * t1527;
t1647 = t1247 * t1405;
t1495 = t1383 * t1567;
t1248 = ((-t1326 * t1380 + t1495) * t1432 + (-t1380 * t1568 - t1383 * t1712) * t1424) * t1400 + ((0.2e1 * t1380 * t1621 + t1383 * t1594) * t1682 + (-t1296 * t1383 - t1380 * t1571) * t1432 - (t1380 * t1350 + t1423 * t1495) * t1621) * t1414 + t1383 * t1503 + (t1380 * t1323 - t1383 * t1450) * t1432 + t1693 * t1383 + t1380 * t1516;
t1646 = t1248 * t1408;
t1493 = t1384 * t1565;
t1249 = ((-t1327 * t1381 + t1493) * t1435 + (-t1381 * t1566 - t1384 * t1711) * t1427) * t1400 + ((0.2e1 * t1381 * t1620 + t1384 * t1592) * t1681 + (-t1297 * t1384 - t1381 * t1570) * t1435 - (t1381 * t1354 + t1426 * t1493) * t1620) * t1414 + t1384 * t1502 + (t1381 * t1324 - t1384 * t1449) * t1435 + t1692 * t1384 + t1381 * t1515;
t1645 = t1249 * t1411;
t1609 = t1420 * t1422;
t1259 = ((-pkin(6) * t1519 - t1325 * t1429) * t1619 + (-t1349 * t1429 - t1430 * t1683) * t1422) * t1414 + t1429 * t1422 * t1468 + (t1685 * t1421 - t1501) * t1609 + ((pkin(6) * t1597 - t1513 * t1421) * t1400 + t1404 * t1679 - t1421 * t1346) * t1431;
t1644 = t1259 * t1405;
t1605 = t1423 * t1425;
t1260 = ((-pkin(6) * t1518 - t1326 * t1432) * t1618 + (-t1353 * t1432 - t1433 * t1682) * t1425) * t1414 + t1432 * t1425 * t1467 + (t1685 * t1424 - t1500) * t1605 + ((pkin(6) * t1595 - t1512 * t1424) * t1400 + t1407 * t1678 - t1350 * t1424) * t1434;
t1643 = t1260 * t1408;
t1601 = t1426 * t1428;
t1261 = ((-pkin(6) * t1517 - t1327 * t1435) * t1617 + (-t1357 * t1435 - t1436 * t1681) * t1428) * t1414 + t1435 * t1428 * t1466 + (t1685 * t1427 - t1499) * t1601 + ((pkin(6) * t1593 - t1511 * t1427) * t1400 + t1410 * t1677 - t1354 * t1427) * t1437;
t1642 = t1261 * t1411;
t1265 = t1716 * t1413 + ((t1416 * t1421 - t1430 * t1707) * pkin(1) + t1710) * t1414;
t1641 = t1265 * t1405;
t1266 = t1715 * t1413 + ((t1416 * t1424 - t1433 * t1706) * pkin(1) + t1709) * t1414;
t1640 = t1266 * t1408;
t1267 = t1714 * t1413 + ((t1416 * t1427 - t1436 * t1705) * pkin(1) + t1708) * t1414;
t1639 = t1267 * t1411;
t1361 = pkin(6) * t1684;
t1362 = pkin(1) * t1668;
t1563 = t1421 * t1676;
t1589 = pkin(5) * t1680;
t1638 = 0.1e1 / ((-t1429 * t1589 + t1361) * t1430 - t1563 * t1684 + t1413 * (-pkin(5) * t1385 + t1362)) * t1405;
t1363 = pkin(1) * t1666;
t1561 = t1424 * t1675;
t1637 = 0.1e1 / ((-t1432 * t1589 + t1361) * t1433 - t1561 * t1684 + t1413 * (-pkin(5) * t1386 + t1363)) * t1408;
t1364 = pkin(1) * t1664;
t1559 = t1427 * t1674;
t1636 = 0.1e1 / ((-t1435 * t1589 + t1361) * t1436 - t1559 * t1684 + t1413 * (-pkin(5) * t1387 + t1364)) * t1411;
t1635 = t1710 * t1405;
t1634 = t1709 * t1408;
t1633 = t1708 * t1411;
t1316 = -t1563 + t1662;
t1632 = t1316 * t1414;
t1317 = -t1561 + t1661;
t1631 = t1317 * t1414;
t1318 = -t1559 + t1660;
t1630 = t1318 * t1414;
t1629 = t1331 * t1405;
t1628 = t1332 * t1408;
t1627 = t1333 * t1411;
t1444 = pkin(2) ^ 2;
t1626 = t1404 * t1444;
t1625 = t1407 * t1444;
t1624 = t1410 * t1444;
t1623 = t1413 * t1414;
t1445 = 0.1e1 / pkin(2);
t1616 = t1413 * t1445;
t1615 = t1414 * t1422;
t1614 = t1414 * t1425;
t1613 = t1414 * t1428;
t1590 = 0.2e1 * pkin(1) * t1416;
t1588 = pkin(6) * t1676;
t1587 = pkin(6) * t1675;
t1586 = pkin(6) * t1674;
t1585 = mrSges(3,1) * t1667;
t1584 = mrSges(3,1) * t1665;
t1583 = mrSges(3,1) * t1663;
t1582 = 0.2e1 * t1623;
t1578 = t1414 * t1668;
t1577 = t1414 * t1666;
t1576 = t1414 * t1664;
t1575 = pkin(2) * t1609;
t1574 = pkin(2) * t1605;
t1573 = pkin(2) * t1601;
t1552 = -0.2e1 * t1588;
t1551 = -0.2e1 * t1587;
t1550 = -0.2e1 * t1586;
t1546 = t1686 * pkin(6);
t1358 = t1388 + pkin(6);
t1347 = pkin(5) + t1668;
t1489 = t1686 * t1422 * t1347;
t1286 = t1358 * t1619 + t1421 * t1489;
t1348 = 0.2e1 * t1385 + pkin(1);
t1289 = t1348 * t1619 + t1489;
t1522 = t1413 * t1615;
t1474 = t1379 * t1522;
t1301 = -t1400 * t1382 + t1382 + t1474;
t1310 = -t1400 * t1385 + t1349;
t1328 = t1388 + t1546;
t1442 = pkin(6) ^ 2;
t1507 = -t1442 + t1626;
t1457 = t1507 * t1431;
t1465 = pkin(6) * t1474;
t1530 = t1347 * t1623;
t1477 = t1382 * t1530;
t1525 = t1421 * t1626;
t1250 = (t1413 * t1379 * t1457 + 0.2e1 * t1301 * t1588) * t1406 + (-t1301 * t1525 + (t1289 * t1379 - t1477) * t1676 + (t1382 * t1310 + t1421 * t1465) * pkin(6)) * t1430 - (t1382 * t1328 + t1465) * t1676 + pkin(6) * (t1286 * t1379 - t1421 * t1477);
t1542 = t1250 * t1638;
t1473 = t1382 * t1522;
t1302 = t1400 * t1379 - t1379 + t1473;
t1464 = pkin(6) * t1473;
t1480 = t1379 * t1530;
t1251 = (-t1507 * t1382 * t1619 + t1302 * t1552) * t1406 + (t1302 * t1525 - (t1289 * t1382 + t1480) * t1676 + pkin(6) * (t1379 * t1310 - t1421 * t1464)) * t1430 - (t1379 * t1328 - t1464) * t1676 - (t1286 * t1382 + t1421 * t1480) * pkin(6);
t1541 = t1251 * t1638;
t1359 = t1389 + pkin(6);
t1351 = pkin(5) + t1666;
t1488 = t1686 * t1425 * t1351;
t1287 = t1359 * t1618 + t1424 * t1488;
t1352 = 0.2e1 * t1386 + pkin(1);
t1290 = t1352 * t1618 + t1488;
t1521 = t1413 * t1614;
t1472 = t1380 * t1521;
t1303 = -t1400 * t1383 + t1383 + t1472;
t1311 = -t1400 * t1386 + t1353;
t1329 = t1389 + t1546;
t1506 = -t1442 + t1625;
t1459 = t1434 * t1506;
t1463 = pkin(6) * t1472;
t1529 = t1351 * t1623;
t1476 = t1383 * t1529;
t1524 = t1424 * t1625;
t1252 = (t1413 * t1380 * t1459 + 0.2e1 * t1303 * t1587) * t1409 + (-t1303 * t1524 + (t1290 * t1380 - t1476) * t1675 + (t1383 * t1311 + t1424 * t1463) * pkin(6)) * t1433 - (t1383 * t1329 + t1463) * t1675 + pkin(6) * (t1287 * t1380 - t1424 * t1476);
t1540 = t1252 * t1637;
t1471 = t1383 * t1521;
t1304 = t1400 * t1380 - t1380 + t1471;
t1462 = pkin(6) * t1471;
t1479 = t1380 * t1529;
t1253 = (-t1506 * t1383 * t1618 + t1304 * t1551) * t1409 + (t1304 * t1524 - (t1290 * t1383 + t1479) * t1675 + pkin(6) * (t1380 * t1311 - t1424 * t1462)) * t1433 - (t1380 * t1329 - t1462) * t1675 - (t1287 * t1383 + t1424 * t1479) * pkin(6);
t1539 = t1253 * t1637;
t1360 = t1390 + pkin(6);
t1355 = pkin(5) + t1664;
t1487 = t1686 * t1428 * t1355;
t1288 = t1360 * t1617 + t1427 * t1487;
t1356 = 0.2e1 * t1387 + pkin(1);
t1291 = t1356 * t1617 + t1487;
t1520 = t1413 * t1613;
t1470 = t1381 * t1520;
t1305 = -t1400 * t1384 + t1384 + t1470;
t1312 = -t1400 * t1387 + t1357;
t1330 = t1390 + t1546;
t1505 = -t1442 + t1624;
t1458 = t1437 * t1505;
t1461 = pkin(6) * t1470;
t1528 = t1355 * t1623;
t1475 = t1384 * t1528;
t1523 = t1427 * t1624;
t1254 = (t1413 * t1381 * t1458 + 0.2e1 * t1305 * t1586) * t1412 + (-t1305 * t1523 + (t1291 * t1381 - t1475) * t1674 + (t1384 * t1312 + t1427 * t1461) * pkin(6)) * t1436 - (t1384 * t1330 + t1461) * t1674 + pkin(6) * (t1288 * t1381 - t1427 * t1475);
t1538 = t1254 * t1636;
t1469 = t1384 * t1520;
t1306 = t1400 * t1381 - t1381 + t1469;
t1460 = pkin(6) * t1469;
t1478 = t1381 * t1528;
t1255 = (-t1505 * t1384 * t1617 + t1306 * t1550) * t1412 + (t1306 * t1523 - (t1291 * t1384 + t1478) * t1674 + pkin(6) * (t1381 * t1312 - t1427 * t1460)) * t1436 - (t1381 * t1330 - t1460) * t1674 - (t1288 * t1384 + t1427 * t1478) * pkin(6);
t1537 = t1255 * t1636;
t1262 = (t1414 * t1431 * t1552 + t1507 * t1422) * t1406 + ((-t1347 * t1619 + t1422 * t1348) * t1676 + t1421 * t1414 * t1457) * t1430 + pkin(6) * (t1422 * t1358 + (-t1347 * t1622 + t1414 * t1676) * t1431);
t1536 = t1262 * t1638;
t1263 = (t1414 * t1434 * t1551 + t1506 * t1425) * t1409 + ((-t1351 * t1618 + t1425 * t1352) * t1675 + t1424 * t1414 * t1459) * t1433 + pkin(6) * (t1425 * t1359 + (-t1351 * t1621 + t1414 * t1675) * t1434);
t1535 = t1263 * t1637;
t1264 = (t1414 * t1437 * t1550 + t1505 * t1428) * t1412 + ((-t1355 * t1617 + t1428 * t1356) * t1674 + t1427 * t1414 * t1458) * t1436 + pkin(6) * (t1428 * t1360 + (-t1355 * t1620 + t1414 * t1674) * t1437);
t1534 = t1264 * t1636;
t1533 = t1445 * t1638;
t1532 = t1445 * t1637;
t1531 = t1445 * t1636;
t1440 = pkin(1) * mrSges(3,1);
t1510 = -t1421 * Ifges(3,5) + t1440;
t1509 = -t1424 * Ifges(3,5) + t1440;
t1508 = -t1427 * Ifges(3,5) + t1440;
t1492 = Ifges(3,3) * t1533;
t1491 = Ifges(3,3) * t1532;
t1490 = Ifges(3,3) * t1531;
t1439 = pkin(1) * mrSges(3,2);
t1274 = ((-mrSges(3,1) * t1673 - t1421 * Ifges(3,6) + t1439) * t1420 + (-mrSges(3,2) * t1673 - t1510) * t1429 - Ifges(3,3) * t1430) * t1413 - ((mrSges(3,1) * t1388 - Ifges(3,5)) * t1420 + (mrSges(3,2) * t1388 - Ifges(3,6)) * t1429) * t1414;
t1486 = t1274 * t1533;
t1275 = ((-mrSges(3,1) * t1672 - t1424 * Ifges(3,6) + t1439) * t1423 + (-mrSges(3,2) * t1672 - t1509) * t1432 - Ifges(3,3) * t1433) * t1413 - ((mrSges(3,1) * t1389 - Ifges(3,5)) * t1423 + (mrSges(3,2) * t1389 - Ifges(3,6)) * t1432) * t1414;
t1485 = t1275 * t1532;
t1276 = ((-mrSges(3,1) * t1671 - t1427 * Ifges(3,6) + t1439) * t1426 + (-mrSges(3,2) * t1671 - t1508) * t1435 - Ifges(3,3) * t1436) * t1413 - ((mrSges(3,1) * t1390 - Ifges(3,5)) * t1426 + (mrSges(3,2) * t1390 - Ifges(3,6)) * t1435) * t1414;
t1484 = t1276 * t1531;
t1483 = t1331 * t1533;
t1482 = t1332 * t1532;
t1481 = t1333 * t1531;
t1456 = t1262 * t1413 * t1533;
t1455 = t1263 * t1413 * t1532;
t1454 = t1264 * t1413 * t1531;
t1453 = Ifges(2,3) - Ifges(2,1) - t1415 + t1698;
t1452 = t1441 * pkin(1) ^ 2 + Ifges(1,3) + t1689 - t1698;
t1438 = pkin(5) * mrSges(3,2);
t1401 = -0.2e1 * t1438;
t1294 = t1318 * t1613 - t1437 * t1701;
t1293 = t1317 * t1614 - t1434 * t1700;
t1292 = t1316 * t1615 - t1431 * t1699;
t1285 = 0.1e1 / (pkin(1) * t1630 + (pkin(5) * t1701 + t1364) * t1413);
t1284 = 0.1e1 / (pkin(1) * t1631 + (pkin(5) * t1700 + t1363) * t1413);
t1283 = 0.1e1 / (pkin(1) * t1632 + (pkin(5) * t1699 + t1362) * t1413);
t1282 = -t1428 * t1701 + (-t1413 * t1664 - t1630) * t1437;
t1281 = -t1425 * t1700 + (-t1413 * t1666 - t1631) * t1434;
t1280 = -t1422 * t1699 + (-t1413 * t1668 - t1632) * t1431;
t1273 = (t1381 * t1318 - t1384 * t1573) * t1413 - t1294 * t1384 - t1381 * t1576;
t1272 = (t1384 * t1318 + t1381 * t1573) * t1413 + t1294 * t1381 - t1384 * t1576;
t1271 = (t1380 * t1317 - t1383 * t1574) * t1413 - t1293 * t1383 - t1380 * t1577;
t1270 = (t1383 * t1317 + t1380 * t1574) * t1413 + t1293 * t1380 - t1383 * t1577;
t1269 = (t1379 * t1316 - t1382 * t1575) * t1413 - t1292 * t1382 - t1379 * t1578;
t1268 = (t1382 * t1316 + t1379 * t1575) * t1413 + t1292 * t1379 - t1382 * t1578;
t1258 = ((-0.2e1 * t1553 + t1704) * t1412 - 0.2e1 * (-Ifges(3,5) * t1435 + t1651) * t1598 + 0.2e1 * t1610 + (t1401 + 0.4e1 * t1654) * t1435 - 0.2e1 * t1583 + t1453) * t1400 + t1714 * t1582 + (t1336 - t1704) * t1412 + 0.2e1 * (-pkin(1) * t1339 + t1651 * t1427 + t1508 * t1435) * t1436 - t1610 + 0.2e1 * (t1438 - t1654) * t1435 + t1427 * t1590 + 0.2e1 * t1583 + t1452;
t1257 = ((-0.2e1 * t1554 + t1703) * t1409 - 0.2e1 * (-Ifges(3,5) * t1432 + t1652) * t1602 + 0.2e1 * t1611 + (t1401 + 0.4e1 * t1655) * t1432 - 0.2e1 * t1584 + t1453) * t1400 + t1715 * t1582 + (t1335 - t1703) * t1409 + 0.2e1 * (-pkin(1) * t1338 + t1652 * t1424 + t1509 * t1432) * t1433 - t1611 + 0.2e1 * (t1438 - t1655) * t1432 + t1424 * t1590 + 0.2e1 * t1584 + t1452;
t1256 = ((-0.2e1 * t1555 + t1702) * t1406 - 0.2e1 * (-Ifges(3,5) * t1429 + t1653) * t1606 + 0.2e1 * t1612 + (t1401 + 0.4e1 * t1656) * t1429 - 0.2e1 * t1585 + t1453) * t1400 + t1716 * t1582 + (t1334 - t1702) * t1406 + 0.2e1 * (-pkin(1) * t1337 + t1653 * t1421 + t1510 * t1429) * t1430 - t1612 + 0.2e1 * (t1438 - t1656) * t1429 + t1421 * t1590 + 0.2e1 * t1585 + t1452;
t1243 = Ifges(3,3) * t1454 + (t1261 * t1627 + t1276 * t1282) * t1285;
t1242 = Ifges(3,3) * t1455 + (t1260 * t1628 + t1275 * t1281) * t1284;
t1241 = Ifges(3,3) * t1456 + (t1259 * t1629 + t1274 * t1280) * t1283;
t1240 = t1333 * t1454 + (t1261 * t1633 + t1267 * t1282) * t1285;
t1239 = t1332 * t1455 + (t1260 * t1634 + t1266 * t1281) * t1284;
t1238 = t1331 * t1456 + (t1259 * t1635 + t1265 * t1280) * t1283;
t1237 = t1276 * t1454 + (t1258 * t1282 + t1261 * t1639) * t1285;
t1236 = t1275 * t1455 + (t1257 * t1281 + t1260 * t1640) * t1284;
t1235 = t1274 * t1456 + (t1256 * t1280 + t1259 * t1641) * t1283;
t1234 = t1255 * t1490 + (t1249 * t1627 + t1273 * t1276) * t1285;
t1233 = t1254 * t1490 + (t1246 * t1627 + t1272 * t1276) * t1285;
t1232 = t1253 * t1491 + (t1248 * t1628 + t1271 * t1275) * t1284;
t1231 = t1252 * t1491 + (t1245 * t1628 + t1270 * t1275) * t1284;
t1230 = t1251 * t1492 + (t1247 * t1629 + t1269 * t1274) * t1283;
t1229 = t1250 * t1492 + (t1244 * t1629 + t1268 * t1274) * t1283;
t1228 = t1255 * t1481 + (t1249 * t1633 + t1267 * t1273) * t1285;
t1227 = t1254 * t1481 + (t1246 * t1633 + t1267 * t1272) * t1285;
t1226 = t1253 * t1482 + (t1248 * t1634 + t1266 * t1271) * t1284;
t1225 = t1252 * t1482 + (t1245 * t1634 + t1266 * t1270) * t1284;
t1224 = t1251 * t1483 + (t1247 * t1635 + t1265 * t1269) * t1283;
t1223 = t1250 * t1483 + (t1244 * t1635 + t1265 * t1268) * t1283;
t1222 = t1255 * t1484 + (t1249 * t1639 + t1258 * t1273) * t1285;
t1221 = t1254 * t1484 + (t1246 * t1639 + t1258 * t1272) * t1285;
t1220 = t1253 * t1485 + (t1248 * t1640 + t1257 * t1271) * t1284;
t1219 = t1252 * t1485 + (t1245 * t1640 + t1257 * t1270) * t1284;
t1218 = t1251 * t1486 + (t1247 * t1641 + t1256 * t1269) * t1283;
t1217 = t1250 * t1486 + (t1244 * t1641 + t1256 * t1268) * t1283;
t1 = [m(4) + (t1222 * t1273 + t1228 * t1645) * t1285 + (t1220 * t1271 + t1226 * t1646) * t1284 + (t1218 * t1269 + t1224 * t1647) * t1283 + (t1230 * t1541 + t1232 * t1539 + t1234 * t1537) * t1445, (t1222 * t1272 + t1228 * t1648) * t1285 + (t1220 * t1270 + t1226 * t1649) * t1284 + (t1218 * t1268 + t1224 * t1650) * t1283 + (t1230 * t1542 + t1232 * t1540 + t1234 * t1538) * t1445, (t1222 * t1282 + t1228 * t1642) * t1285 + (t1220 * t1281 + t1226 * t1643) * t1284 + (t1218 * t1280 + t1224 * t1644) * t1283 + (t1230 * t1536 + t1232 * t1535 + t1234 * t1534) * t1616; (t1221 * t1273 + t1227 * t1645) * t1285 + (t1219 * t1271 + t1225 * t1646) * t1284 + (t1217 * t1269 + t1223 * t1647) * t1283 + (t1229 * t1541 + t1231 * t1539 + t1233 * t1537) * t1445, m(4) + (t1221 * t1272 + t1227 * t1648) * t1285 + (t1219 * t1270 + t1225 * t1649) * t1284 + (t1217 * t1268 + t1223 * t1650) * t1283 + (t1229 * t1542 + t1231 * t1540 + t1233 * t1538) * t1445, (t1221 * t1282 + t1227 * t1642) * t1285 + (t1219 * t1281 + t1225 * t1643) * t1284 + (t1217 * t1280 + t1223 * t1644) * t1283 + (t1229 * t1536 + t1231 * t1535 + t1233 * t1534) * t1616; (t1237 * t1273 + t1240 * t1645) * t1285 + (t1236 * t1271 + t1239 * t1646) * t1284 + (t1235 * t1269 + t1238 * t1647) * t1283 + (t1241 * t1541 + t1242 * t1539 + t1243 * t1537) * t1445, (t1237 * t1272 + t1240 * t1648) * t1285 + (t1236 * t1270 + t1239 * t1649) * t1284 + (t1235 * t1268 + t1238 * t1650) * t1283 + (t1241 * t1542 + t1242 * t1540 + t1243 * t1538) * t1445, m(4) + (t1237 * t1282 + t1240 * t1642) * t1285 + (t1236 * t1281 + t1239 * t1643) * t1284 + (t1235 * t1280 + t1238 * t1644) * t1283 + (t1241 * t1536 + t1242 * t1535 + t1243 * t1534) * t1616;];
MX  = t1;
