% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRP1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*(3+1)/2x13]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:31
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRP1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A0_inertia_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:12:24
% EndTime: 2019-05-03 15:12:52
% DurationCPUTime: 30.76s
% Computational Cost: add. (75131->918), mult. (166297->1453), div. (1272->9), fcn. (65150->26), ass. (0->661)
t1601 = legFrame(3,3);
t1551 = cos(t1601);
t1864 = qJ(3,3) * t1551;
t1737 = pkin(2) * t1864;
t1512 = -0.2e1 * t1737;
t1548 = sin(t1601);
t1617 = qJ(3,3) ^ 2;
t1626 = pkin(2) ^ 2;
t1555 = -t1617 + t1626;
t1433 = t1548 * t1555 + t1512;
t1521 = t1548 * qJ(3,3);
t1515 = pkin(2) * t1521;
t1611 = cos(qJ(1,3));
t1599 = t1626 / 0.2e1;
t1717 = t1599 - t1617 / 0.2e1;
t1605 = sin(qJ(1,3));
t1901 = -0.2e1 * t1605;
t1379 = t1433 * t1611 + (t1551 * t1717 + t1515) * t1901;
t1509 = 0.2e1 * t1515;
t1385 = (t1551 * t1555 + t1509) * t1611 + t1433 * t1605;
t1524 = t1548 * pkin(2);
t1451 = t1524 - t1864;
t1862 = qJ(3,3) * t1605;
t1886 = pkin(2) * t1611;
t1482 = t1862 - t1886;
t1592 = qJ(1,3) + qJ(2,3);
t1545 = cos(t1592);
t1536 = t1545 ^ 2;
t1542 = sin(t1592);
t1627 = (pkin(1) ^ 2);
t1889 = pkin(2) * t1551;
t1449 = t1521 - t1889;
t1450 = t1524 + t1864;
t1589 = t1611 ^ 2;
t1776 = t1605 * t1611;
t1638 = -t1449 * t1776 + t1450 * t1589;
t1639 = -t1449 * t1589 - t1450 * t1776;
t1785 = t1545 * t1551;
t1597 = 1 + t1627;
t1790 = t1542 * t1597;
t1791 = t1542 * t1545;
t1331 = (t1548 * t1790 - t1785) * qJ(3,3) + (t1379 * t1536 + t1385 * t1791 + t1451 * t1482) * pkin(1) + ((-t1638 + t1524) * t1545 + t1639 * t1542) * t1627;
t1539 = 0.1e1 + t1555;
t1430 = t1539 * t1548 + t1512;
t1794 = t1539 * t1551;
t1382 = (t1509 + t1794) * t1611 + t1605 * t1430;
t1427 = t1430 * t1611;
t1904 = -0.2e1 * t1515;
t1346 = (t1542 * t1548 - t1785) * qJ(3,3) + (-t1427 * t1536 - t1382 * t1791 + (t1548 * t1626 + t1548 - t1737) * t1611 + (-(t1904 - t1794) * t1536 - qJ(3,3) * t1451) * t1605) * pkin(1);
t1753 = t1331 + t1346;
t1602 = legFrame(2,3);
t1552 = cos(t1602);
t1871 = qJ(3,2) * t1552;
t1738 = pkin(2) * t1871;
t1513 = -0.2e1 * t1738;
t1549 = sin(t1602);
t1618 = qJ(3,2) ^ 2;
t1557 = -t1618 + t1626;
t1434 = t1549 * t1557 + t1513;
t1522 = t1549 * qJ(3,2);
t1516 = pkin(2) * t1522;
t1613 = cos(qJ(1,2));
t1716 = t1599 - t1618 / 0.2e1;
t1607 = sin(qJ(1,2));
t1900 = -0.2e1 * t1607;
t1380 = t1434 * t1613 + (t1552 * t1716 + t1516) * t1900;
t1510 = 0.2e1 * t1516;
t1386 = (t1552 * t1557 + t1510) * t1613 + t1434 * t1607;
t1525 = t1549 * pkin(2);
t1455 = t1525 - t1871;
t1869 = qJ(3,2) * t1607;
t1885 = pkin(2) * t1613;
t1483 = t1869 - t1885;
t1593 = qJ(1,2) + qJ(2,2);
t1546 = cos(t1593);
t1537 = t1546 ^ 2;
t1543 = sin(t1593);
t1888 = pkin(2) * t1552;
t1453 = t1522 - t1888;
t1454 = t1525 + t1871;
t1590 = t1613 ^ 2;
t1775 = t1607 * t1613;
t1636 = -t1453 * t1775 + t1454 * t1590;
t1637 = -t1453 * t1590 - t1454 * t1775;
t1784 = t1546 * t1552;
t1788 = t1543 * t1597;
t1789 = t1543 * t1546;
t1332 = (t1549 * t1788 - t1784) * qJ(3,2) + (t1380 * t1537 + t1386 * t1789 + t1455 * t1483) * pkin(1) + ((-t1636 + t1525) * t1546 + t1637 * t1543) * t1627;
t1540 = 0.1e1 + t1557;
t1431 = t1540 * t1549 + t1513;
t1793 = t1540 * t1552;
t1383 = (t1510 + t1793) * t1613 + t1607 * t1431;
t1428 = t1431 * t1613;
t1903 = -0.2e1 * t1516;
t1347 = (t1543 * t1549 - t1784) * qJ(3,2) + (-t1428 * t1537 - t1383 * t1789 + (t1549 * t1626 + t1549 - t1738) * t1613 + (-(t1903 - t1793) * t1537 - qJ(3,2) * t1455) * t1607) * pkin(1);
t1752 = t1332 + t1347;
t1603 = legFrame(1,3);
t1553 = cos(t1603);
t1878 = qJ(3,1) * t1553;
t1741 = pkin(2) * t1878;
t1514 = -0.2e1 * t1741;
t1550 = sin(t1603);
t1619 = qJ(3,1) ^ 2;
t1559 = -t1619 + t1626;
t1435 = t1550 * t1559 + t1514;
t1523 = t1550 * qJ(3,1);
t1517 = pkin(2) * t1523;
t1615 = cos(qJ(1,1));
t1715 = t1599 - t1619 / 0.2e1;
t1609 = sin(qJ(1,1));
t1899 = -0.2e1 * t1609;
t1381 = t1435 * t1615 + (t1553 * t1715 + t1517) * t1899;
t1511 = 0.2e1 * t1517;
t1387 = (t1553 * t1559 + t1511) * t1615 + t1435 * t1609;
t1526 = t1550 * pkin(2);
t1459 = t1526 - t1878;
t1876 = qJ(3,1) * t1609;
t1884 = pkin(2) * t1615;
t1484 = t1876 - t1884;
t1594 = qJ(1,1) + qJ(2,1);
t1547 = cos(t1594);
t1538 = t1547 ^ 2;
t1544 = sin(t1594);
t1887 = pkin(2) * t1553;
t1457 = t1523 - t1887;
t1458 = t1526 + t1878;
t1591 = t1615 ^ 2;
t1774 = t1609 * t1615;
t1634 = -t1457 * t1774 + t1458 * t1591;
t1635 = -t1457 * t1591 - t1458 * t1774;
t1783 = t1547 * t1553;
t1786 = t1544 * t1597;
t1787 = t1544 * t1547;
t1333 = (t1550 * t1786 - t1783) * qJ(3,1) + (t1381 * t1538 + t1387 * t1787 + t1459 * t1484) * pkin(1) + ((-t1634 + t1526) * t1547 + t1635 * t1544) * t1627;
t1541 = 0.1e1 + t1559;
t1432 = t1541 * t1550 + t1514;
t1792 = t1541 * t1553;
t1384 = (t1511 + t1792) * t1615 + t1609 * t1432;
t1429 = t1432 * t1615;
t1902 = -0.2e1 * t1517;
t1348 = (t1544 * t1550 - t1783) * qJ(3,1) + (-t1429 * t1538 - t1384 * t1787 + (t1550 * t1626 + t1550 - t1741) * t1615 + (-(t1902 - t1792) * t1538 - qJ(3,1) * t1459) * t1609) * pkin(1);
t1751 = t1333 + t1348;
t1448 = t1521 + t1889;
t1711 = t1542 * t1864;
t1797 = t1521 * t1545;
t1334 = -t1597 * t1711 - t1797 + (t1379 * t1791 - t1385 * t1536 - t1448 * t1482) * pkin(1) + ((t1639 - t1889) * t1545 + t1638 * t1542) * t1627;
t1518 = 0.1e1 / 0.2e1 + t1717;
t1800 = t1448 * t1605;
t1343 = -t1711 - t1797 + (t1382 * t1536 - (t1427 + (t1518 * t1551 + t1515) * t1901) * t1791 - (t1551 * t1626 + t1515 + t1551) * t1611 + qJ(3,3) * t1800) * pkin(1);
t1750 = t1334 + t1343;
t1452 = t1522 + t1888;
t1710 = t1543 * t1871;
t1796 = t1522 * t1546;
t1335 = -t1597 * t1710 - t1796 + (t1380 * t1789 - t1386 * t1537 - t1452 * t1483) * pkin(1) + ((t1637 - t1888) * t1546 + t1636 * t1543) * t1627;
t1519 = 0.1e1 / 0.2e1 + t1716;
t1799 = t1452 * t1607;
t1344 = -t1710 - t1796 + (t1383 * t1537 - (t1428 + (t1519 * t1552 + t1516) * t1900) * t1789 - (t1552 * t1626 + t1516 + t1552) * t1613 + qJ(3,2) * t1799) * pkin(1);
t1749 = t1335 + t1344;
t1456 = t1523 + t1887;
t1709 = t1544 * t1878;
t1795 = t1523 * t1547;
t1336 = -t1597 * t1709 - t1795 + (t1381 * t1787 - t1387 * t1538 - t1456 * t1484) * pkin(1) + ((t1635 - t1887) * t1547 + t1634 * t1544) * t1627;
t1520 = 0.1e1 / 0.2e1 + t1715;
t1798 = t1456 * t1609;
t1345 = -t1709 - t1795 + (t1384 * t1538 - (t1429 + (t1520 * t1553 + t1517) * t1899) * t1787 - (t1553 * t1626 + t1517 + t1553) * t1615 + qJ(3,1) * t1798) * pkin(1);
t1748 = t1336 + t1345;
t1896 = pkin(2) * qJ(3,3);
t1685 = t1776 * t1896;
t1631 = -t1539 * t1589 + 0.2e1 * t1685;
t1764 = -0.1e1 / 0.2e1 - t1626 / 0.2e1;
t1596 = 2 + t1627;
t1865 = qJ(3,3) * t1542;
t1916 = 0.2e1 * pkin(1);
t1909 = -(-(pkin(1) * (t1539 * t1776 + (0.2e1 * t1589 - 0.1e1) * t1896) * t1542 + t1862) * t1545 + t1611 * t1865) * t1916 - t1617 * t1596;
t1358 = (0.2e1 * (t1617 / 0.2e1 - t1631 + t1764) * t1536 + t1631) * t1627 + t1909;
t1352 = 0.1e1 / t1358;
t1554 = t1617 + t1626;
t1925 = t1352 * t1554;
t1897 = pkin(2) * qJ(3,2);
t1686 = t1775 * t1897;
t1632 = -t1540 * t1590 + 0.2e1 * t1686;
t1872 = qJ(3,2) * t1543;
t1908 = -(-(pkin(1) * (t1540 * t1775 + (0.2e1 * t1590 - 0.1e1) * t1897) * t1543 + t1869) * t1546 + t1613 * t1872) * t1916 - t1618 * t1596;
t1359 = (0.2e1 * (t1618 / 0.2e1 - t1632 + t1764) * t1537 + t1632) * t1627 + t1908;
t1354 = 0.1e1 / t1359;
t1556 = t1618 + t1626;
t1924 = t1354 * t1556;
t1898 = pkin(2) * qJ(3,1);
t1687 = t1774 * t1898;
t1633 = -t1541 * t1591 + 0.2e1 * t1687;
t1879 = qJ(3,1) * t1544;
t1907 = -(-(pkin(1) * (t1541 * t1774 + (0.2e1 * t1591 - 0.1e1) * t1898) * t1544 + t1876) * t1547 + t1615 * t1879) * t1916 - t1619 * t1596;
t1360 = (0.2e1 * (t1619 / 0.2e1 - t1633 + t1764) * t1538 + t1633) * t1627 + t1907;
t1356 = 0.1e1 / t1360;
t1558 = t1619 + t1626;
t1923 = t1356 * t1558;
t1616 = xP(3);
t1566 = sin(t1616);
t1567 = cos(t1616);
t1621 = koppelP(2,2);
t1624 = koppelP(2,1);
t1443 = t1566 * t1624 + t1567 * t1621;
t1818 = t1354 * t1443;
t1283 = t1332 * t1818;
t1446 = -t1566 * t1621 + t1567 * t1624;
t1817 = t1354 * t1446;
t1285 = t1335 * t1817;
t1247 = -t1283 + t1285;
t1314 = t1344 * t1817;
t1317 = t1347 * t1818;
t1281 = -t1317 + t1314;
t1922 = t1247 + t1281;
t1622 = koppelP(1,2);
t1625 = koppelP(1,1);
t1444 = t1566 * t1625 + t1567 * t1622;
t1812 = t1356 * t1444;
t1284 = t1333 * t1812;
t1447 = -t1566 * t1622 + t1567 * t1625;
t1811 = t1356 * t1447;
t1286 = t1336 * t1811;
t1248 = -t1284 + t1286;
t1315 = t1345 * t1811;
t1318 = t1348 * t1812;
t1282 = -t1318 + t1315;
t1921 = t1248 + t1282;
t1620 = koppelP(3,2);
t1623 = koppelP(3,1);
t1442 = t1566 * t1623 + t1567 * t1620;
t1824 = t1352 * t1442;
t1287 = t1331 * t1824;
t1445 = -t1566 * t1620 + t1567 * t1623;
t1823 = t1352 * t1445;
t1288 = t1334 * t1823;
t1249 = -t1287 + t1288;
t1313 = t1343 * t1823;
t1316 = t1346 * t1824;
t1280 = -t1316 + t1313;
t1920 = t1249 + t1280;
t1905 = 0.2e1 * pkin(2);
t1915 = t1922 * t1905;
t1914 = t1921 * t1905;
t1913 = t1920 * t1905;
t1560 = 0.2e1 * t1617 + t1627;
t1779 = t1589 * t1627;
t1912 = t1560 - t1779;
t1561 = 0.2e1 * t1618 + t1627;
t1778 = t1590 * t1627;
t1911 = t1561 - t1778;
t1562 = 0.2e1 * t1619 + t1627;
t1777 = t1591 * t1627;
t1910 = t1562 - t1777;
t1604 = sin(qJ(2,3));
t1895 = pkin(1) * t1604;
t1606 = sin(qJ(2,2));
t1894 = pkin(1) * t1606;
t1608 = sin(qJ(2,1));
t1893 = pkin(1) * t1608;
t1610 = cos(qJ(2,3));
t1892 = pkin(1) * t1610;
t1612 = cos(qJ(2,2));
t1891 = pkin(1) * t1612;
t1614 = cos(qJ(2,1));
t1890 = pkin(1) * t1614;
t1883 = pkin(2) * t1627;
t1882 = qJ(3,1) * t1921;
t1261 = t1748 * t1356;
t1881 = qJ(3,1) * t1261;
t1273 = t1751 * t1356;
t1880 = qJ(3,1) * t1273;
t1877 = qJ(3,1) * t1608;
t1575 = qJ(3,1) * t1622;
t1576 = qJ(3,1) * t1625;
t1875 = qJ(3,2) * t1922;
t1260 = t1749 * t1354;
t1874 = qJ(3,2) * t1260;
t1272 = t1752 * t1354;
t1873 = qJ(3,2) * t1272;
t1870 = qJ(3,2) * t1606;
t1572 = qJ(3,2) * t1621;
t1573 = qJ(3,2) * t1624;
t1868 = qJ(3,3) * t1920;
t1259 = t1750 * t1352;
t1867 = qJ(3,3) * t1259;
t1271 = t1753 * t1352;
t1866 = qJ(3,3) * t1271;
t1863 = qJ(3,3) * t1604;
t1569 = qJ(3,3) * t1620;
t1570 = qJ(3,3) * t1623;
t1736 = pkin(2) * t1569;
t1528 = -0.2e1 * t1736;
t1568 = t1617 * t1623;
t1770 = t1623 * t1626;
t1684 = -t1568 + t1770;
t1466 = t1528 + t1684;
t1734 = pkin(2) * t1570;
t1531 = 0.2e1 * t1734;
t1577 = t1626 * t1620;
t1467 = -t1617 * t1620 + t1531 + t1577;
t1400 = t1466 * t1567 - t1467 * t1566;
t1401 = t1466 * t1566 + t1467 * t1567;
t1364 = -t1400 * t1548 + t1401 * t1551;
t1583 = pkin(2) * t1623;
t1491 = t1569 + t1583;
t1580 = pkin(2) * t1620;
t1496 = -t1570 + t1580;
t1415 = t1491 * t1567 - t1496 * t1566;
t1418 = t1491 * t1566 + t1496 * t1567;
t1367 = t1415 * t1551 + t1418 * t1548;
t1492 = -t1569 + t1583;
t1495 = t1570 + t1580;
t1416 = t1492 * t1567 - t1495 * t1566;
t1417 = t1492 * t1566 + t1495 * t1567;
t1368 = t1416 * t1551 + t1417 * t1548;
t1370 = -t1415 * t1548 + t1418 * t1551;
t1397 = t1442 * t1548 + t1445 * t1551;
t1769 = t1623 * t1627;
t1485 = pkin(2) * t1769 - t1569;
t1486 = t1620 * t1883 + t1570;
t1660 = t1400 * t1551 + t1548 * t1401;
t1232 = ((t1485 * t1567 - t1486 * t1566) * t1551 + (t1485 * t1566 + t1486 * t1567) * t1548) * t1545 + qJ(3,3) * t1397 * t1790 + ((-t1367 * t1589 - t1370 * t1776) * t1545 + (-t1367 * t1776 + t1370 * t1589) * t1542) * t1627 + (-(t1605 * t1364 - t1611 * t1660) * t1536 + (t1364 * t1611 + t1605 * t1660) * t1791 + t1368 * t1482) * pkin(1);
t1763 = 0.1e1 / 0.4e1 + t1626 / 0.4e1;
t1349 = 0.1e1 / ((0.4e1 * (-t1518 * t1589 + t1685 - t1617 / 0.4e1 + t1763) * t1536 - t1631) * t1627 - t1909);
t1860 = t1232 * t1349;
t1735 = pkin(2) * t1572;
t1530 = -0.2e1 * t1735;
t1571 = t1618 * t1624;
t1768 = t1624 * t1626;
t1683 = -t1571 + t1768;
t1468 = t1530 + t1683;
t1733 = pkin(2) * t1573;
t1532 = 0.2e1 * t1733;
t1578 = t1626 * t1621;
t1469 = -t1618 * t1621 + t1532 + t1578;
t1402 = t1468 * t1567 - t1469 * t1566;
t1403 = t1468 * t1566 + t1469 * t1567;
t1365 = -t1402 * t1549 + t1403 * t1552;
t1584 = pkin(2) * t1624;
t1497 = t1572 + t1584;
t1581 = pkin(2) * t1621;
t1502 = -t1573 + t1581;
t1419 = t1497 * t1567 - t1502 * t1566;
t1422 = t1497 * t1566 + t1502 * t1567;
t1371 = t1419 * t1552 + t1422 * t1549;
t1498 = -t1572 + t1584;
t1501 = t1573 + t1581;
t1420 = t1498 * t1567 - t1501 * t1566;
t1421 = t1498 * t1566 + t1501 * t1567;
t1372 = t1420 * t1552 + t1421 * t1549;
t1374 = -t1419 * t1549 + t1422 * t1552;
t1398 = t1443 * t1549 + t1446 * t1552;
t1767 = t1624 * t1627;
t1487 = pkin(2) * t1767 - t1572;
t1488 = t1621 * t1883 + t1573;
t1659 = t1402 * t1552 + t1549 * t1403;
t1233 = ((t1487 * t1567 - t1488 * t1566) * t1552 + (t1487 * t1566 + t1488 * t1567) * t1549) * t1546 + qJ(3,2) * t1398 * t1788 + ((-t1371 * t1590 - t1374 * t1775) * t1546 + (-t1371 * t1775 + t1374 * t1590) * t1543) * t1627 + (-(t1607 * t1365 - t1613 * t1659) * t1537 + (t1365 * t1613 + t1607 * t1659) * t1789 + t1372 * t1483) * pkin(1);
t1350 = 0.1e1 / ((0.4e1 * (-t1519 * t1590 + t1686 - t1618 / 0.4e1 + t1763) * t1537 - t1632) * t1627 - t1908);
t1859 = t1233 * t1350;
t1740 = pkin(2) * t1575;
t1534 = -0.2e1 * t1740;
t1574 = t1619 * t1625;
t1766 = t1625 * t1626;
t1682 = -t1574 + t1766;
t1470 = t1534 + t1682;
t1739 = pkin(2) * t1576;
t1535 = 0.2e1 * t1739;
t1579 = t1626 * t1622;
t1471 = -t1619 * t1622 + t1535 + t1579;
t1404 = t1470 * t1567 - t1471 * t1566;
t1405 = t1470 * t1566 + t1471 * t1567;
t1366 = -t1404 * t1550 + t1405 * t1553;
t1585 = pkin(2) * t1625;
t1503 = t1575 + t1585;
t1582 = pkin(2) * t1622;
t1508 = -t1576 + t1582;
t1423 = t1503 * t1567 - t1508 * t1566;
t1426 = t1503 * t1566 + t1508 * t1567;
t1375 = t1423 * t1553 + t1426 * t1550;
t1504 = -t1575 + t1585;
t1507 = t1576 + t1582;
t1424 = t1504 * t1567 - t1507 * t1566;
t1425 = t1504 * t1566 + t1507 * t1567;
t1376 = t1424 * t1553 + t1425 * t1550;
t1378 = -t1423 * t1550 + t1426 * t1553;
t1399 = t1444 * t1550 + t1447 * t1553;
t1765 = t1625 * t1627;
t1489 = pkin(2) * t1765 - t1575;
t1490 = t1622 * t1883 + t1576;
t1658 = t1404 * t1553 + t1550 * t1405;
t1234 = ((t1489 * t1567 - t1490 * t1566) * t1553 + (t1489 * t1566 + t1490 * t1567) * t1550) * t1547 + qJ(3,1) * t1399 * t1786 + ((-t1375 * t1591 - t1378 * t1774) * t1547 + (-t1375 * t1774 + t1378 * t1591) * t1544) * t1627 + (-(t1609 * t1366 - t1615 * t1658) * t1538 + (t1366 * t1615 + t1609 * t1658) * t1787 + t1376 * t1484) * pkin(1);
t1351 = 0.1e1 / ((0.4e1 * (-t1520 * t1591 + t1687 - t1619 / 0.4e1 + t1763) * t1538 - t1633) * t1627 - t1907);
t1858 = t1234 * t1351;
t1595 = 0.1e1 + t1626;
t1436 = t1531 + (t1595 - t1617) * t1620;
t1437 = t1595 * t1623 + t1528 - t1568;
t1388 = t1436 * t1567 + t1437 * t1566;
t1782 = t1566 * t1436;
t1389 = t1437 * t1567 - t1782;
t1361 = t1388 * t1551 - t1389 * t1548;
t1394 = t1442 * t1551 - t1445 * t1548;
t1475 = t1623 - t1736 + t1770;
t1476 = t1577 + t1620 + t1734;
t1527 = 0.2e1 * t1736;
t1803 = t1388 * t1548;
t1806 = t1368 * t1605;
t1250 = (-t1394 * t1545 + t1397 * t1542) * qJ(3,3) + (((((t1527 - t1623 - t1684) * t1567 + t1782) * t1551 - t1803) * t1611 + t1361 * t1605) * t1536 - (t1361 * t1611 + (t1389 * t1551 + t1803) * t1605) * t1791 - ((-t1475 * t1567 + t1476 * t1566) * t1551 - (t1475 * t1566 + t1476 * t1567) * t1548) * t1611 - qJ(3,3) * t1806) * pkin(1);
t1857 = t1250 * t1349;
t1438 = t1532 + (t1595 - t1618) * t1621;
t1439 = t1595 * t1624 + t1530 - t1571;
t1390 = t1438 * t1567 + t1439 * t1566;
t1781 = t1566 * t1438;
t1391 = t1439 * t1567 - t1781;
t1362 = t1390 * t1552 - t1391 * t1549;
t1395 = t1443 * t1552 - t1446 * t1549;
t1477 = t1624 - t1735 + t1768;
t1478 = t1578 + t1621 + t1733;
t1529 = 0.2e1 * t1735;
t1802 = t1390 * t1549;
t1805 = t1372 * t1607;
t1251 = (-t1395 * t1546 + t1398 * t1543) * qJ(3,2) + (((((t1529 - t1624 - t1683) * t1567 + t1781) * t1552 - t1802) * t1613 + t1362 * t1607) * t1537 - (t1362 * t1613 + (t1391 * t1552 + t1802) * t1607) * t1789 - ((-t1477 * t1567 + t1478 * t1566) * t1552 - (t1477 * t1566 + t1478 * t1567) * t1549) * t1613 - qJ(3,2) * t1805) * pkin(1);
t1856 = t1251 * t1350;
t1440 = t1535 + (t1595 - t1619) * t1622;
t1441 = t1595 * t1625 + t1534 - t1574;
t1392 = t1440 * t1567 + t1441 * t1566;
t1780 = t1566 * t1440;
t1393 = t1441 * t1567 - t1780;
t1363 = t1392 * t1553 - t1393 * t1550;
t1396 = t1444 * t1553 - t1447 * t1550;
t1479 = t1625 - t1740 + t1766;
t1480 = t1579 + t1622 + t1739;
t1533 = 0.2e1 * t1740;
t1801 = t1392 * t1550;
t1804 = t1376 * t1609;
t1252 = (-t1396 * t1547 + t1399 * t1544) * qJ(3,1) + (((((t1533 - t1625 - t1682) * t1567 + t1780) * t1553 - t1801) * t1615 + t1363 * t1609) * t1538 - (t1363 * t1615 + (t1393 * t1553 + t1801) * t1609) * t1787 - ((-t1479 * t1567 + t1480 * t1566) * t1553 - (t1479 * t1566 + t1480 * t1567) * t1550) * t1615 - qJ(3,1) * t1804) * pkin(1);
t1855 = t1252 * t1351;
t1854 = t1334 * t1352;
t1853 = t1335 * t1354;
t1852 = t1336 * t1356;
t1406 = t1448 * t1611 + t1451 * t1605;
t1407 = t1451 * t1611 - t1800;
t1705 = t1627 * t1776;
t1472 = t1548 * t1705;
t1678 = t1551 * t1705;
t1681 = -0.2e1 * t1617 - t1779;
t1337 = (t1551 * t1912 + t1472 + t1904) * t1545 + (t1548 * t1681 + t1512 - t1678) * t1542 + (-t1407 * t1791 + t1406 * t1536 - 0.2e1 * t1521 * t1611 + (-t1862 - t1886) * t1551) * pkin(1);
t1851 = t1337 * t1352;
t1408 = t1452 * t1613 + t1455 * t1607;
t1409 = t1455 * t1613 - t1799;
t1704 = t1627 * t1775;
t1473 = t1549 * t1704;
t1677 = t1552 * t1704;
t1680 = -0.2e1 * t1618 - t1778;
t1338 = (t1552 * t1911 + t1473 + t1903) * t1546 + (t1549 * t1680 + t1513 - t1677) * t1543 + (-t1409 * t1789 + t1408 * t1537 - 0.2e1 * t1522 * t1613 + (-t1869 - t1885) * t1552) * pkin(1);
t1850 = t1338 * t1354;
t1410 = t1456 * t1615 + t1459 * t1609;
t1411 = t1459 * t1615 - t1798;
t1703 = t1627 * t1774;
t1474 = t1550 * t1703;
t1676 = t1553 * t1703;
t1679 = -0.2e1 * t1619 - t1777;
t1339 = (t1553 * t1910 + t1474 + t1902) * t1547 + (t1550 * t1679 + t1514 - t1676) * t1544 + (-t1411 * t1787 + t1410 * t1538 - 0.2e1 * t1523 * t1615 + (-t1876 - t1884) * t1553) * pkin(1);
t1849 = t1339 * t1356;
t1340 = (-t1548 * t1912 + t1512 + t1678) * t1545 + (t1551 * t1681 + t1472 + t1509) * t1542 + (-t1407 * t1536 - t1406 * t1791 + (t1524 - 0.2e1 * t1864) * t1611 + t1605 * t1521) * pkin(1);
t1848 = t1340 * t1352;
t1341 = (-t1549 * t1911 + t1513 + t1677) * t1546 + (t1552 * t1680 + t1473 + t1510) * t1543 + (-t1409 * t1537 - t1408 * t1789 + (t1525 - 0.2e1 * t1871) * t1613 + t1607 * t1522) * pkin(1);
t1847 = t1341 * t1354;
t1342 = (-t1550 * t1910 + t1514 + t1676) * t1547 + (t1553 * t1679 + t1474 + t1511) * t1544 + (-t1411 * t1538 - t1410 * t1787 + (t1526 - 0.2e1 * t1878) * t1615 + t1609 * t1523) * pkin(1);
t1846 = t1342 * t1356;
t1845 = t1343 * t1352;
t1844 = t1343 * t1610;
t1843 = t1344 * t1354;
t1842 = t1344 * t1612;
t1841 = t1345 * t1356;
t1840 = t1345 * t1614;
t1839 = t1346 * t1352;
t1353 = 0.1e1 / t1358 ^ 2;
t1838 = t1346 * t1353;
t1837 = t1346 * t1610;
t1836 = t1347 * t1354;
t1355 = 0.1e1 / t1359 ^ 2;
t1835 = t1347 * t1355;
t1834 = t1347 * t1612;
t1833 = t1348 * t1356;
t1357 = 0.1e1 / t1360 ^ 2;
t1832 = t1348 * t1357;
t1831 = t1348 * t1614;
t1830 = t1349 * t1604;
t1829 = t1349 * t1610;
t1828 = t1350 * t1606;
t1827 = t1350 * t1612;
t1826 = t1351 * t1608;
t1825 = t1351 * t1614;
t1819 = t1352 * t1627;
t1813 = t1354 * t1627;
t1807 = t1356 * t1627;
t1726 = 0.2e1 * t1839;
t1744 = t1352 * t1905;
t1759 = pkin(2) * t1726 + t1331 * t1744;
t1725 = 0.2e1 * t1836;
t1743 = t1354 * t1905;
t1758 = pkin(2) * t1725 + t1332 * t1743;
t1724 = 0.2e1 * t1833;
t1742 = t1356 * t1905;
t1757 = pkin(2) * t1724 + t1333 * t1742;
t1729 = 0.2e1 * t1845;
t1756 = pkin(2) * t1729 + t1334 * t1744;
t1728 = 0.2e1 * t1843;
t1755 = pkin(2) * t1728 + t1335 * t1743;
t1727 = 0.2e1 * t1841;
t1754 = pkin(2) * t1727 + t1336 * t1742;
t1307 = t1339 * t1811;
t1310 = t1342 * t1812;
t1253 = -t1310 + t1307;
t1308 = t1337 * t1823;
t1311 = t1340 * t1824;
t1254 = -t1311 + t1308;
t1309 = t1338 * t1817;
t1312 = t1341 * t1818;
t1255 = -t1312 + t1309;
t1747 = 0.2e1 * t1892;
t1746 = 0.2e1 * t1891;
t1745 = 0.2e1 * t1890;
t1732 = 0.2e1 * t1857;
t1731 = 0.2e1 * t1856;
t1730 = 0.2e1 * t1855;
t1723 = t1343 + t1334 / 0.2e1;
t1722 = t1344 + t1335 / 0.2e1;
t1721 = t1345 + t1336 / 0.2e1;
t1720 = t1346 + t1331 / 0.2e1;
t1719 = t1347 + t1332 / 0.2e1;
t1718 = t1348 + t1333 / 0.2e1;
t1714 = t1352 * t1895;
t1713 = t1354 * t1894;
t1712 = t1356 * t1893;
t1708 = t1352 * t1857;
t1707 = t1354 * t1856;
t1706 = t1356 * t1855;
t1702 = t1877 * t1916;
t1701 = t1870 * t1916;
t1700 = t1863 * t1916;
t1699 = t1753 * t1925;
t1698 = t1750 * t1925;
t1697 = t1752 * t1924;
t1696 = t1749 * t1924;
t1695 = t1751 * t1923;
t1694 = t1748 * t1923;
t1693 = t1352 * t1747;
t1692 = t1354 * t1746;
t1691 = t1356 * t1745;
t1690 = t1351 * (t1252 + t1234);
t1689 = (t1250 + t1232) * t1349;
t1688 = (t1251 + t1233) * t1350;
t1675 = pkin(1) * t1844 - t1337;
t1674 = pkin(1) * t1842 - t1338;
t1673 = pkin(1) * t1840 - t1339;
t1672 = pkin(1) * t1837 - t1340;
t1671 = pkin(1) * t1834 - t1341;
t1670 = pkin(1) * t1831 - t1342;
t1669 = (-0.2e1 * t1317 + 0.2e1 * t1314 + t1247) * t1251 + t1233 * t1281;
t1668 = (-0.2e1 * t1318 + 0.2e1 * t1315 + t1248) * t1252 + t1234 * t1282;
t1667 = (-0.2e1 * t1316 + 0.2e1 * t1313 + t1249) * t1250 + t1232 * t1280;
t1262 = (t1334 + 0.2e1 * t1343) * t1352;
t1666 = t1262 * t1352 + t1334 * t1353;
t1263 = (t1335 + 0.2e1 * t1344) * t1354;
t1665 = t1263 * t1354 + t1335 * t1355;
t1264 = (t1336 + 0.2e1 * t1345) * t1356;
t1664 = t1264 * t1356 + t1336 * t1357;
t1268 = (t1331 + 0.2e1 * t1346) * t1352;
t1663 = t1268 * t1352 + t1331 * t1353;
t1269 = (t1332 + 0.2e1 * t1347) * t1354;
t1662 = t1269 * t1354 + t1332 * t1355;
t1270 = (t1333 + 0.2e1 * t1348) * t1356;
t1661 = t1270 * t1356 + t1333 * t1357;
t1657 = -t1280 * t1892 + t1254;
t1656 = -t1281 * t1891 + t1255;
t1655 = -t1282 * t1890 + t1253;
t1654 = -pkin(2) * t1254 + t1554 * t1920;
t1653 = -pkin(2) * t1255 + t1556 * t1922;
t1652 = -pkin(2) * t1253 + t1558 * t1921;
t1651 = pkin(1) * (pkin(2) * t1614 + t1877);
t1650 = pkin(1) * (pkin(2) * t1612 + t1870);
t1649 = pkin(1) * (pkin(2) * t1610 + t1863);
t1645 = t1232 * t1845 + t1250 * t1262;
t1644 = t1232 * t1839 + t1250 * t1268;
t1643 = t1233 * t1843 + t1251 * t1263;
t1642 = t1233 * t1836 + t1251 * t1269;
t1641 = t1234 * t1841 + t1252 * t1264;
t1640 = t1234 * t1833 + t1252 * t1270;
t1630 = t1268 * t1845 + t1334 * t1838;
t1629 = t1269 * t1843 + t1335 * t1835;
t1628 = t1270 * t1841 + t1336 * t1832;
t1506 = 0.2e1 * t1576 + t1582;
t1505 = -0.2e1 * t1575 + t1585;
t1500 = 0.2e1 * t1573 + t1581;
t1499 = -0.2e1 * t1572 + t1584;
t1494 = 0.2e1 * t1570 + t1580;
t1493 = -0.2e1 * t1569 + t1583;
t1481 = t1566 ^ 2 + t1567 ^ 2;
t1465 = t1562 * t1622 - 0.2e1 * t1739;
t1464 = t1533 + 0.2e1 * t1574 + t1765;
t1463 = t1561 * t1621 - 0.2e1 * t1733;
t1462 = t1529 + 0.2e1 * t1571 + t1767;
t1461 = t1560 * t1620 - 0.2e1 * t1734;
t1460 = t1527 + 0.2e1 * t1568 + t1769;
t1377 = -t1424 * t1550 + t1425 * t1553;
t1373 = -t1420 * t1549 + t1421 * t1552;
t1369 = -t1416 * t1548 + t1417 * t1551;
t1267 = t1718 * t1356;
t1266 = t1719 * t1354;
t1265 = t1720 * t1352;
t1258 = t1721 * t1356;
t1257 = t1722 * t1354;
t1256 = t1723 * t1352;
t1237 = ((-t1464 * t1567 + t1465 * t1566) * t1553 - t1550 * (t1464 * t1566 + t1465 * t1567)) * t1547 + 0.2e1 * t1376 * t1879 + ((t1396 * t1774 + t1399 * t1591) * t1547 + (-t1396 * t1591 + t1399 * t1774) * t1544) * t1627 + ((-t1376 * t1615 + t1609 * t1377) * t1538 - (t1377 * t1615 + t1804) * t1787 + ((t1505 * t1567 - t1506 * t1566) * t1553 + t1550 * (t1505 * t1566 + t1506 * t1567)) * t1615 + t1399 * t1876) * pkin(1);
t1236 = ((-t1462 * t1567 + t1463 * t1566) * t1552 - t1549 * (t1462 * t1566 + t1463 * t1567)) * t1546 + 0.2e1 * t1372 * t1872 + ((t1395 * t1775 + t1398 * t1590) * t1546 + (-t1395 * t1590 + t1398 * t1775) * t1543) * t1627 + ((-t1372 * t1613 + t1607 * t1373) * t1537 - (t1373 * t1613 + t1805) * t1789 + ((t1499 * t1567 - t1500 * t1566) * t1552 + t1549 * (t1499 * t1566 + t1500 * t1567)) * t1613 + t1398 * t1869) * pkin(1);
t1235 = ((-t1460 * t1567 + t1461 * t1566) * t1551 - t1548 * (t1460 * t1566 + t1461 * t1567)) * t1545 + 0.2e1 * t1368 * t1865 + ((t1394 * t1776 + t1397 * t1589) * t1545 + (-t1394 * t1589 + t1397 * t1776) * t1542) * t1627 + ((-t1368 * t1611 + t1605 * t1369) * t1536 - (t1369 * t1611 + t1806) * t1791 + ((t1493 * t1567 - t1494 * t1566) * t1551 + (t1493 * t1566 + t1494 * t1567) * t1548) * t1611 + t1397 * t1862) * pkin(1);
t1231 = t1348 * t1712 + 0.2e1 * t1880;
t1230 = t1347 * t1713 + 0.2e1 * t1873;
t1229 = t1346 * t1714 + 0.2e1 * t1866;
t1228 = t1345 * t1712 + 0.2e1 * t1881;
t1227 = t1344 * t1713 + 0.2e1 * t1874;
t1226 = t1343 * t1714 + 0.2e1 * t1867;
t1225 = (-pkin(2) * t1751 - t1670) * t1356;
t1224 = t1356 * t1670 + t1757;
t1223 = (-pkin(2) * t1752 - t1671) * t1354;
t1222 = t1354 * t1671 + t1758;
t1221 = (-pkin(2) * t1753 - t1672) * t1352;
t1220 = t1352 * t1672 + t1759;
t1219 = (-pkin(2) * t1748 - t1673) * t1356;
t1218 = t1356 * t1673 + t1754;
t1217 = (-pkin(2) * t1749 - t1674) * t1354;
t1216 = t1354 * t1674 + t1755;
t1215 = (-pkin(2) * t1750 - t1675) * t1352;
t1214 = t1352 * t1675 + t1756;
t1213 = t1266 * t1894 + t1873;
t1212 = t1267 * t1893 + t1880;
t1211 = t1265 * t1895 + t1866;
t1210 = t1257 * t1894 + t1874;
t1209 = t1258 * t1893 + t1881;
t1208 = t1256 * t1895 + t1867;
t1201 = -t1284 / 0.2e1 + t1286 / 0.2e1 + t1282;
t1200 = -t1283 / 0.2e1 + t1285 / 0.2e1 + t1281;
t1199 = -t1287 / 0.2e1 + t1288 / 0.2e1 + t1280;
t1198 = t1267 * t1745 + t1757 - t1846;
t1197 = t1266 * t1746 + t1758 - t1847;
t1196 = t1265 * t1747 + t1759 - t1848;
t1195 = t1258 * t1745 + t1754 - t1849;
t1194 = t1257 * t1746 + t1755 - t1850;
t1193 = t1256 * t1747 + t1756 - t1851;
t1192 = (-pkin(2) * t1342 + t1348 * t1651) * t1356 + t1695;
t1191 = (-pkin(2) * t1341 + t1347 * t1650) * t1354 + t1697;
t1190 = (-pkin(2) * t1340 + t1346 * t1649) * t1352 + t1699;
t1189 = (-pkin(2) * t1339 + t1345 * t1651) * t1356 + t1694;
t1188 = (-pkin(2) * t1338 + t1344 * t1650) * t1354 + t1696;
t1187 = (-pkin(2) * t1337 + t1343 * t1649) * t1352 + t1698;
t1186 = (-t1342 / 0.2e1 + t1718 * pkin(2)) * t1691 + t1267 * t1702 + t1348 * t1807 - pkin(2) * t1846 + t1695;
t1185 = (-t1341 / 0.2e1 + t1719 * pkin(2)) * t1692 + t1266 * t1701 + t1347 * t1813 - pkin(2) * t1847 + t1697;
t1184 = (-t1340 / 0.2e1 + t1720 * pkin(2)) * t1693 + t1265 * t1700 + t1346 * t1819 - pkin(2) * t1848 + t1699;
t1183 = (-t1339 / 0.2e1 + t1721 * pkin(2)) * t1691 + t1258 * t1702 + t1345 * t1807 - pkin(2) * t1849 + t1694;
t1182 = (-t1338 / 0.2e1 + t1722 * pkin(2)) * t1692 + t1257 * t1701 + t1344 * t1813 - pkin(2) * t1850 + t1696;
t1181 = (-t1337 / 0.2e1 + t1723 * pkin(2)) * t1693 + t1256 * t1700 + t1343 * t1819 - pkin(2) * t1851 + t1698;
t1 = [t1346 ^ 2 * t1353 + t1347 ^ 2 * t1355 + t1348 ^ 2 * t1357, 0, 0, t1271 ^ 2 + t1272 ^ 2 + t1273 ^ 2, (t1661 * t1831 + t1662 * t1834 + t1663 * t1837) * pkin(1), (-t1346 * t1604 * t1663 - t1347 * t1606 * t1662 - t1348 * t1608 * t1661) * pkin(1), (t1198 * t1348 + t1224 * t1333 - t1273 * t1342) * t1356 + (t1197 * t1347 + t1222 * t1332 - t1272 * t1341) * t1354 + (t1196 * t1346 + t1220 * t1331 - t1271 * t1340) * t1352, t1331 * t1352 * t1229 + t1332 * t1354 * t1230 + t1333 * t1356 * t1231 + t1211 * t1726 + t1212 * t1724 + t1213 * t1725, (t1186 * t1348 + t1192 * t1333 + t1225 * t1342) * t1356 + (t1185 * t1347 + t1191 * t1332 + t1223 * t1341) * t1354 + (t1184 * t1346 + t1190 * t1331 + t1221 * t1340) * t1352, 0, 0, 0, t1481; t1343 * t1838 + t1344 * t1835 + t1345 * t1832, 0, 0, t1259 * t1271 + t1260 * t1272 + t1261 * t1273, (t1610 * t1630 + t1612 * t1629 + t1614 * t1628) * pkin(1), (-t1604 * t1630 - t1606 * t1629 - t1608 * t1628) * pkin(1), (t1198 * t1345 + t1224 * t1336 - t1273 * t1339) * t1356 + (t1197 * t1344 + t1222 * t1335 - t1272 * t1338) * t1354 + (t1196 * t1343 + t1220 * t1334 - t1271 * t1337) * t1352, t1211 * t1729 + t1212 * t1727 + t1213 * t1728 + t1229 * t1854 + t1230 * t1853 + t1231 * t1852, (t1186 * t1345 + t1192 * t1336 + t1225 * t1339) * t1356 + (t1185 * t1344 + t1191 * t1335 + t1223 * t1338) * t1354 + (t1184 * t1343 + t1190 * t1334 + t1221 * t1337) * t1352, 0, 0, 0, 0; t1343 ^ 2 * t1353 + t1344 ^ 2 * t1355 + t1345 ^ 2 * t1357, 0, 0, t1259 ^ 2 + t1260 ^ 2 + t1261 ^ 2, (t1664 * t1840 + t1665 * t1842 + t1666 * t1844) * pkin(1), (-t1343 * t1604 * t1666 - t1344 * t1606 * t1665 - t1345 * t1608 * t1664) * pkin(1), (t1195 * t1345 + t1218 * t1336 - t1261 * t1339) * t1356 + (t1194 * t1344 + t1216 * t1335 - t1260 * t1338) * t1354 + (t1193 * t1343 + t1214 * t1334 - t1259 * t1337) * t1352, t1208 * t1729 + t1209 * t1727 + t1210 * t1728 + t1226 * t1854 + t1227 * t1853 + t1228 * t1852, (t1183 * t1345 + t1189 * t1336 + t1219 * t1339) * t1356 + (t1182 * t1344 + t1188 * t1335 + t1217 * t1338) * t1354 + (t1181 * t1343 + t1187 * t1334 + t1215 * t1337) * t1352, 0, 0, 0, t1481; t1346 * t1708 + t1347 * t1707 + t1348 * t1706, 0, 0, t1271 * t1689 + t1272 * t1688 + t1273 * t1690, (t1640 * t1825 + t1642 * t1827 + t1644 * t1829) * pkin(1), (-t1640 * t1826 - t1642 * t1828 - t1644 * t1830) * pkin(1), (t1198 * t1252 + t1224 * t1234 - t1237 * t1273) * t1351 + (t1197 * t1251 + t1222 * t1233 - t1236 * t1272) * t1350 + (t1196 * t1250 + t1220 * t1232 - t1235 * t1271) * t1349, t1211 * t1732 + t1212 * t1730 + t1213 * t1731 + t1229 * t1860 + t1230 * t1859 + t1231 * t1858, (t1186 * t1252 + t1192 * t1234 + t1225 * t1237) * t1351 + (t1185 * t1251 + t1191 * t1233 + t1223 * t1236) * t1350 + (t1184 * t1250 + t1190 * t1232 + t1221 * t1235) * t1349, 0, -t1566, -t1567, 0; t1343 * t1708 + t1344 * t1707 + t1345 * t1706, 0, 0, t1259 * t1689 + t1260 * t1688 + t1261 * t1690, (t1641 * t1825 + t1643 * t1827 + t1645 * t1829) * pkin(1), (-t1641 * t1826 - t1643 * t1828 - t1645 * t1830) * pkin(1), (t1195 * t1252 + t1218 * t1234 - t1237 * t1261) * t1351 + (t1194 * t1251 + t1216 * t1233 - t1236 * t1260) * t1350 + (t1193 * t1250 + t1214 * t1232 - t1235 * t1259) * t1349, t1208 * t1732 + t1209 * t1730 + t1210 * t1731 + t1226 * t1860 + t1227 * t1859 + t1228 * t1858, (t1183 * t1252 + t1189 * t1234 + t1219 * t1237) * t1351 + (t1182 * t1251 + t1188 * t1233 + t1217 * t1236) * t1350 + (t1181 * t1250 + t1187 * t1232 + t1215 * t1235) * t1349, 0, t1567, -t1566, 0; t1280 * t1857 + t1281 * t1856 + t1282 * t1855, 0, 0, t1688 * t1922 + t1689 * t1920 + t1690 * t1921, (t1667 * t1829 + t1668 * t1825 + t1669 * t1827) * pkin(1), (-t1667 * t1830 - t1668 * t1826 - t1669 * t1828) * pkin(1), (t1252 * (t1201 * t1745 - t1253 + t1914) + t1234 * (-t1655 + t1914) - t1237 * t1921) * t1351 + (t1251 * (t1200 * t1746 - t1255 + t1915) + t1233 * (-t1656 + t1915) - t1236 * t1922) * t1350 + (t1250 * (t1199 * t1747 - t1254 + t1913) + t1232 * (-t1657 + t1913) - t1235 * t1920) * t1349, (t1282 * t1893 + 0.2e1 * t1882) * t1858 + (t1281 * t1894 + 0.2e1 * t1875) * t1859 + (t1280 * t1895 + 0.2e1 * t1868) * t1860 + (t1201 * t1893 + t1882) * t1730 + (t1200 * t1894 + t1875) * t1731 + (t1199 * t1895 + t1868) * t1732, ((t1310 / 0.2e1 - t1307 / 0.2e1 + (t1282 + t1248 / 0.2e1) * pkin(2)) * t1745 + t1201 * t1702 + t1627 * t1282 + t1652) * t1855 + (t1282 * t1651 + t1652) * t1858 + t1237 * t1351 * (-pkin(2) * t1921 + t1655) + ((t1312 / 0.2e1 - t1309 / 0.2e1 + (t1281 + t1247 / 0.2e1) * pkin(2)) * t1746 + t1200 * t1701 + t1627 * t1281 + t1653) * t1856 + (t1281 * t1650 + t1653) * t1859 + t1236 * t1350 * (-pkin(2) * t1922 + t1656) + ((t1311 / 0.2e1 - t1308 / 0.2e1 + (t1280 + t1249 / 0.2e1) * pkin(2)) * t1747 + t1199 * t1700 + t1627 * t1280 + t1654) * t1857 + (t1280 * t1649 + t1654) * t1860 + t1235 * t1349 * (-pkin(2) * t1920 + t1657), 1, 0, 0, 0;];
tau_reg  = t1;