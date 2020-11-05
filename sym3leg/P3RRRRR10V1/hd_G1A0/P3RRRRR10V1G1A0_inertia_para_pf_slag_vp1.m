% Calculate inertia matrix for parallel robot
% P3RRRRR10V1G1A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 22:09
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:36:09
% EndTime: 2020-08-06 21:36:20
% DurationCPUTime: 10.16s
% Computational Cost: add. (13467->550), mult. (30177->886), div. (648->13), fcn. (20088->32), ass. (0->392)
t1767 = m(3) * rSges(3,1);
t1663 = (pkin(5) * t1767);
t1815 = -2 * t1663;
t1513 = cos(qJ(3,1));
t1766 = m(3) * rSges(3,2);
t1446 = rSges(3,1) * t1766 - Icges(3,4);
t1504 = sin(qJ(3,1));
t1672 = t1504 * t1446;
t1526 = rSges(3,2) ^ 2;
t1528 = rSges(3,1) ^ 2;
t1462 = -t1526 + t1528;
t1405 = t1462 * m(3) - Icges(3,1) + Icges(3,2);
t1490 = t1513 ^ 2;
t1803 = t1405 * t1490;
t1814 = -0.2e1 * t1513 * t1672 + t1803;
t1510 = cos(qJ(3,2));
t1501 = sin(qJ(3,2));
t1675 = t1501 * t1446;
t1487 = t1510 ^ 2;
t1804 = t1405 * t1487;
t1813 = -0.2e1 * t1510 * t1675 + t1804;
t1507 = cos(qJ(3,3));
t1498 = sin(qJ(3,3));
t1678 = t1498 * t1446;
t1484 = t1507 ^ 2;
t1805 = t1405 * t1484;
t1812 = -0.2e1 * t1507 * t1678 + t1805;
t1508 = cos(qJ(2,3));
t1486 = t1508 ^ 2;
t1525 = rSges(3,3) ^ 2;
t1527 = rSges(2,2) ^ 2;
t1529 = rSges(2,1) ^ 2;
t1738 = Icges(2,1) + Icges(3,2);
t1802 = Icges(2,2) + Icges(3,3) - t1738 + (-t1527 + t1529) * m(2) + (-t1525 + t1526) * m(3);
t1811 = -(t1802 + t1812) * t1486 + t1498 * t1815;
t1511 = cos(qJ(2,2));
t1489 = t1511 ^ 2;
t1810 = -(t1802 + t1813) * t1489 + t1501 * t1815;
t1514 = cos(qJ(2,1));
t1492 = t1514 ^ 2;
t1809 = -(t1802 + t1814) * t1492 + t1504 * t1815;
t1494 = cos(pkin(3));
t1483 = t1494 ^ 2;
t1745 = t1483 * pkin(6);
t1493 = sin(pkin(3));
t1689 = t1493 * t1494;
t1765 = rSges(3,3) * m(3);
t1444 = rSges(3,2) * t1765 - Icges(3,6);
t1445 = rSges(3,1) * t1765 - Icges(3,5);
t1375 = -t1444 * t1507 - t1445 * t1498;
t1376 = -t1444 * t1510 - t1445 * t1501;
t1377 = -t1444 * t1513 - t1445 * t1504;
t1505 = sin(qJ(2,1));
t1670 = t1505 * t1514;
t1603 = t1504 * t1670;
t1752 = pkin(2) * t1513;
t1651 = t1493 * t1752;
t1664 = t1514 * t1490;
t1686 = t1493 * t1504;
t1470 = t1505 * pkin(6);
t1791 = t1470 + pkin(1);
t1801 = (pkin(2) * t1664 + t1513 * t1791) * t1494 - t1603 * t1651 + pkin(6) * (t1514 - 0.1e1) * (t1514 + 0.1e1) * t1686;
t1502 = sin(qJ(2,2));
t1673 = t1502 * t1511;
t1605 = t1501 * t1673;
t1753 = pkin(2) * t1510;
t1652 = t1493 * t1753;
t1666 = t1511 * t1487;
t1687 = t1493 * t1501;
t1469 = t1502 * pkin(6);
t1792 = t1469 + pkin(1);
t1800 = (pkin(2) * t1666 + t1510 * t1792) * t1494 - t1605 * t1652 + pkin(6) * (t1511 - 0.1e1) * (t1511 + 0.1e1) * t1687;
t1499 = sin(qJ(2,3));
t1676 = t1499 * t1508;
t1607 = t1498 * t1676;
t1754 = pkin(2) * t1507;
t1653 = t1493 * t1754;
t1668 = t1508 * t1484;
t1688 = t1493 * t1498;
t1468 = t1499 * pkin(6);
t1793 = t1468 + pkin(1);
t1799 = (pkin(2) * t1668 + t1507 * t1793) * t1494 - t1607 * t1653 + pkin(6) * (t1508 - 0.1e1) * (t1508 + 0.1e1) * t1688;
t1516 = rSges(2,3) + pkin(5);
t1762 = m(2) * t1516;
t1542 = -rSges(2,1) * t1762 + Icges(2,5) + t1446;
t1578 = -rSges(2,2) * t1762 + pkin(5) * t1765 + Icges(2,6);
t1657 = t1498 * t1766;
t1769 = -0.2e1 * t1446;
t1798 = (t1578 - t1375) * t1508 + ((-t1405 * t1498 - t1663) * t1507 + t1484 * t1769 + pkin(5) * t1657 + t1542) * t1499;
t1656 = t1501 * t1766;
t1797 = (t1578 - t1376) * t1511 + ((-t1405 * t1501 - t1663) * t1510 + t1487 * t1769 + pkin(5) * t1656 + t1542) * t1502;
t1655 = t1504 * t1766;
t1796 = (t1578 - t1377) * t1514 + ((-t1405 * t1504 - t1663) * t1513 + t1490 * t1769 + pkin(5) * t1655 + t1542) * t1505;
t1795 = -0.2e1 * pkin(6);
t1794 = 0.2e1 * pkin(6);
t1742 = t1504 * pkin(2);
t1451 = pkin(1) * t1742;
t1665 = t1513 * t1514;
t1639 = t1505 * t1752;
t1395 = -t1514 * pkin(6) + t1639;
t1703 = t1395 * t1494;
t1356 = 0.1e1 / (pkin(1) * t1703 + (-t1451 + (pkin(2) * t1665 + t1470) * pkin(5)) * t1493);
t1743 = t1501 * pkin(2);
t1450 = pkin(1) * t1743;
t1667 = t1510 * t1511;
t1641 = t1502 * t1753;
t1394 = -t1511 * pkin(6) + t1641;
t1704 = t1394 * t1494;
t1355 = 0.1e1 / (pkin(1) * t1704 + (-t1450 + (pkin(2) * t1667 + t1469) * pkin(5)) * t1493);
t1744 = t1498 * pkin(2);
t1449 = pkin(1) * t1744;
t1669 = t1507 * t1508;
t1643 = t1499 * t1754;
t1393 = -t1508 * pkin(6) + t1643;
t1705 = t1393 * t1494;
t1354 = 0.1e1 / (pkin(1) * t1705 + (-t1449 + (pkin(2) * t1669 + t1468) * pkin(5)) * t1493);
t1760 = pkin(1) * t1499;
t1759 = pkin(1) * t1502;
t1758 = pkin(1) * t1505;
t1768 = m(2) * rSges(2,2);
t1736 = rSges(2,1) * t1768 - Icges(2,4);
t1495 = legFrame(3,3);
t1452 = sin(t1495);
t1455 = cos(t1495);
t1500 = sin(qJ(1,3));
t1509 = cos(qJ(1,3));
t1381 = t1452 * t1500 - t1509 * t1455;
t1530 = pkin(6) ^ 2;
t1532 = pkin(2) ^ 2;
t1591 = t1484 * t1532 - t1530;
t1790 = t1381 * t1591;
t1496 = legFrame(2,3);
t1453 = sin(t1496);
t1456 = cos(t1496);
t1503 = sin(qJ(1,2));
t1512 = cos(qJ(1,2));
t1382 = t1453 * t1503 - t1512 * t1456;
t1590 = t1487 * t1532 - t1530;
t1789 = t1382 * t1590;
t1497 = legFrame(1,3);
t1454 = sin(t1497);
t1457 = cos(t1497);
t1506 = sin(qJ(1,1));
t1515 = cos(qJ(1,1));
t1383 = t1454 * t1506 - t1515 * t1457;
t1589 = t1490 * t1532 - t1530;
t1788 = t1383 * t1589;
t1401 = (t1486 - 0.2e1) * t1744 - pkin(5);
t1425 = t1498 * pkin(5) + pkin(2);
t1757 = pkin(2) * t1484;
t1598 = t1425 - 0.2e1 * t1757;
t1774 = (-pkin(6) * t1607 - t1401 * t1507) * t1689 + t1669 * t1745 + (t1598 * t1483 - t1425 + t1757) * t1499;
t1402 = (t1489 - 0.2e1) * t1743 - pkin(5);
t1428 = t1501 * pkin(5) + pkin(2);
t1756 = pkin(2) * t1487;
t1597 = t1428 - 0.2e1 * t1756;
t1773 = (-pkin(6) * t1605 - t1402 * t1510) * t1689 + t1667 * t1745 + (t1597 * t1483 - t1428 + t1756) * t1502;
t1403 = (t1492 - 0.2e1) * t1742 - pkin(5);
t1431 = t1504 * pkin(5) + pkin(2);
t1755 = pkin(2) * t1490;
t1596 = t1431 - 0.2e1 * t1755;
t1772 = (-pkin(6) * t1603 - t1403 * t1513) * t1689 + t1665 * t1745 + (t1596 * t1483 - t1431 + t1755) * t1505;
t1771 = -0.2e1 * pkin(1);
t1770 = (t1483 - 0.1e1) * t1794;
t1519 = pkin(1) * rSges(3,1);
t1764 = pkin(1) * rSges(3,2);
t1763 = t1405 / 0.2e1;
t1761 = pkin(1) * t1494;
t1751 = pkin(5) * t1508;
t1750 = pkin(5) * t1511;
t1749 = pkin(5) * t1514;
t1747 = (t1527 + t1529) * m(2);
t1737 = Icges(3,1) + Icges(2,3);
t1384 = t1452 * t1509 + t1500 * t1455;
t1617 = t1384 * t1688;
t1315 = t1774 * t1381 + t1799 * t1384 - t1617 * t1760;
t1485 = 0.1e1 / t1507;
t1735 = t1315 * t1485;
t1385 = t1453 * t1512 + t1503 * t1456;
t1616 = t1385 * t1687;
t1316 = t1773 * t1382 + t1800 * t1385 - t1616 * t1759;
t1488 = 0.1e1 / t1510;
t1734 = t1316 * t1488;
t1386 = t1454 * t1515 + t1506 * t1457;
t1615 = t1386 * t1686;
t1317 = t1772 * t1383 + t1801 * t1386 - t1615 * t1758;
t1491 = 0.1e1 / t1513;
t1733 = t1317 * t1491;
t1620 = t1381 * t1688;
t1318 = -t1799 * t1381 + t1774 * t1384 + t1620 * t1760;
t1732 = t1318 * t1485;
t1619 = t1382 * t1687;
t1319 = -t1800 * t1382 + t1773 * t1385 + t1619 * t1759;
t1731 = t1319 * t1488;
t1618 = t1383 * t1686;
t1320 = -t1801 * t1383 + t1772 * t1386 + t1618 * t1758;
t1730 = t1320 * t1491;
t1520 = m(2) * rSges(2,1);
t1415 = -t1520 + t1657;
t1424 = -t1765 + t1768;
t1582 = (t1525 + t1526) * m(3) + t1747 + t1737;
t1327 = t1798 * t1493 + ((-t1424 * t1499 + (t1507 * t1767 - t1415) * t1508) * pkin(1) + t1582 + t1812) * t1494;
t1729 = t1327 * t1485;
t1416 = -t1520 + t1656;
t1328 = t1797 * t1493 + ((-t1424 * t1502 + (t1510 * t1767 - t1416) * t1511) * pkin(1) + t1582 + t1813) * t1494;
t1728 = t1328 * t1488;
t1417 = -t1520 + t1655;
t1329 = t1796 * t1493 + ((-t1424 * t1505 + (t1513 * t1767 - t1417) * t1514) * pkin(1) + t1582 + t1814) * t1494;
t1727 = t1329 * t1491;
t1426 = pkin(5) + t1744;
t1434 = pkin(6) + t1760;
t1685 = t1493 * t1508;
t1610 = t1494 * t1685;
t1614 = (t1494 + 0.1e1) * (t1494 - 0.1e1) * t1532;
t1681 = t1494 * t1499;
t1701 = t1426 * t1493;
t1339 = -t1499 * t1614 * t1668 + (t1426 * t1610 + t1486 * t1770 + t1434 - t1745) * t1754 - ((-t1483 * t1468 + t1793) * t1508 - t1681 * t1701) * pkin(6);
t1533 = 0.1e1 / pkin(2);
t1726 = t1339 * t1533;
t1429 = pkin(5) + t1743;
t1435 = pkin(6) + t1759;
t1684 = t1493 * t1511;
t1609 = t1494 * t1684;
t1680 = t1494 * t1502;
t1700 = t1429 * t1493;
t1340 = -t1502 * t1614 * t1666 + (t1429 * t1609 + t1489 * t1770 + t1435 - t1745) * t1753 - ((-t1483 * t1469 + t1792) * t1511 - t1680 * t1700) * pkin(6);
t1725 = t1340 * t1533;
t1432 = pkin(5) + t1742;
t1436 = pkin(6) + t1758;
t1683 = t1493 * t1514;
t1608 = t1494 * t1683;
t1679 = t1494 * t1505;
t1699 = t1432 * t1493;
t1341 = -t1505 * t1614 * t1664 + (t1432 * t1608 + t1492 * t1770 + t1436 - t1745) * t1752 - ((-t1483 * t1470 + t1791) * t1514 - t1679 * t1699) * pkin(6);
t1724 = t1341 * t1533;
t1521 = 0.2e1 * qJ(3,3);
t1660 = t1526 + t1528;
t1574 = Icges(2,3) + t1747 + (0.2e1 * t1525 + t1660) * m(3) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(3,1) / 0.2e1;
t1345 = cos(t1521) * t1763 - t1446 * sin(t1521) + t1574;
t1723 = t1345 * t1485;
t1522 = 0.2e1 * qJ(3,2);
t1346 = cos(t1522) * t1763 - t1446 * sin(t1522) + t1574;
t1722 = t1346 * t1488;
t1523 = 0.2e1 * qJ(3,1);
t1347 = cos(t1523) * t1763 - t1446 * sin(t1523) + t1574;
t1721 = t1347 * t1491;
t1447 = pkin(6) * t1761;
t1692 = t1493 * (-pkin(5) * t1468 + t1449);
t1720 = 0.1e1 / ((-pkin(5) * t1653 + t1447) * t1508 - t1643 * t1761 + t1692) * t1485;
t1691 = t1493 * (-pkin(5) * t1469 + t1450);
t1719 = 0.1e1 / ((-pkin(5) * t1652 + t1447) * t1511 - t1641 * t1761 + t1691) * t1488;
t1690 = t1493 * (-pkin(5) * t1470 + t1451);
t1718 = 0.1e1 / ((-pkin(5) * t1651 + t1447) * t1514 - t1639 * t1761 + t1690) * t1491;
t1717 = 0.1e1 / ((pkin(1) * t1681 + pkin(5) * t1685) * t1754 - t1508 * t1447 - t1692) * t1485;
t1716 = 0.1e1 / ((pkin(1) * t1680 + pkin(5) * t1684) * t1753 - t1447 * t1511 - t1691) * t1488;
t1715 = 0.1e1 / ((pkin(1) * t1679 + pkin(5) * t1683) * t1752 - t1514 * t1447 - t1690) * t1491;
t1714 = t1354 * (t1393 * t1493 + t1494 * t1744);
t1713 = t1355 * (t1394 * t1493 + t1494 * t1743);
t1712 = t1356 * (t1395 * t1493 + t1494 * t1742);
t1708 = t1375 * t1485;
t1707 = t1376 * t1488;
t1706 = t1377 * t1491;
t1418 = t1660 * m(3) + Icges(3,3);
t1702 = t1418 * t1533;
t1698 = t1444 * t1498;
t1697 = t1444 * t1501;
t1696 = t1444 * t1504;
t1695 = t1445 * t1499;
t1694 = t1445 * t1502;
t1693 = t1445 * t1505;
t1682 = t1493 * t1533;
t1677 = t1498 * t1508;
t1674 = t1501 * t1511;
t1671 = t1504 * t1514;
t1662 = pkin(5) * t1766;
t1661 = 0.2e1 * t1736;
t1531 = pkin(5) ^ 2;
t1659 = pkin(1) ^ 2 + t1531;
t1658 = t1424 * t1771;
t1654 = -0.2e1 * t1689;
t1650 = t1494 * t1754;
t1649 = t1494 * t1753;
t1648 = t1494 * t1752;
t1613 = t1509 * t1701;
t1360 = t1500 * t1434 - t1499 * t1613;
t1427 = 0.2e1 * t1468 + pkin(1);
t1369 = t1500 * t1427 - t1613;
t1606 = t1500 * t1701;
t1548 = t1434 * t1509 + t1499 * t1606;
t1551 = t1427 * t1509 + t1606;
t1573 = t1591 * t1384;
t1585 = t1384 * t1650;
t1321 = (t1585 * t1795 + t1790) * t1486 + ((t1452 * t1369 - t1455 * t1551) * t1754 + t1573 * t1681) * t1508 + (t1452 * t1360 - t1455 * t1548 + t1585) * pkin(6);
t1632 = t1321 * t1720;
t1612 = t1512 * t1700;
t1361 = t1503 * t1435 - t1502 * t1612;
t1430 = 0.2e1 * t1469 + pkin(1);
t1370 = t1503 * t1430 - t1612;
t1604 = t1503 * t1700;
t1547 = t1435 * t1512 + t1502 * t1604;
t1550 = t1430 * t1512 + t1604;
t1572 = t1590 * t1385;
t1584 = t1385 * t1649;
t1322 = (t1584 * t1795 + t1789) * t1489 + ((t1453 * t1370 - t1456 * t1550) * t1753 + t1572 * t1680) * t1511 + (t1453 * t1361 - t1456 * t1547 + t1584) * pkin(6);
t1631 = t1322 * t1719;
t1611 = t1515 * t1699;
t1362 = t1506 * t1436 - t1505 * t1611;
t1433 = 0.2e1 * t1470 + pkin(1);
t1371 = t1506 * t1433 - t1611;
t1602 = t1506 * t1699;
t1546 = t1436 * t1515 + t1505 * t1602;
t1549 = t1433 * t1515 + t1602;
t1571 = t1589 * t1386;
t1583 = t1386 * t1648;
t1323 = (t1583 * t1795 + t1788) * t1492 + ((t1454 * t1371 - t1457 * t1549) * t1752 + t1571 * t1679) * t1514 + (t1454 * t1362 - t1457 * t1546 + t1583) * pkin(6);
t1630 = t1323 * t1718;
t1588 = t1381 * t1650;
t1324 = (t1588 * t1794 + t1573) * t1486 + ((t1369 * t1455 + t1452 * t1551) * t1754 - t1681 * t1790) * t1508 + (t1360 * t1455 + t1452 * t1548 - t1588) * pkin(6);
t1629 = t1324 * t1720;
t1587 = t1382 * t1649;
t1325 = (t1587 * t1794 + t1572) * t1489 + ((t1370 * t1456 + t1453 * t1550) * t1753 - t1680 * t1789) * t1511 + (t1361 * t1456 + t1453 * t1547 - t1587) * pkin(6);
t1628 = t1325 * t1719;
t1586 = t1383 * t1648;
t1326 = (t1586 * t1794 + t1571) * t1492 + ((t1371 * t1457 + t1454 * t1549) * t1752 - t1679 * t1788) * t1514 + (t1362 * t1457 + t1454 * t1546 - t1586) * pkin(6);
t1627 = t1326 * t1718;
t1626 = t1339 * t1717;
t1625 = t1340 * t1716;
t1624 = t1341 * t1715;
t1342 = (pkin(6) * t1610 + t1401 * t1483 + pkin(5) + (-t1486 + 0.1e1) * t1744) * t1507 - t1793 * t1677 + (t1598 * t1689 + t1677 * t1745) * t1499;
t1623 = t1342 * t1717;
t1343 = (pkin(6) * t1609 + t1402 * t1483 + pkin(5) + (-t1489 + 0.1e1) * t1743) * t1510 - t1792 * t1674 + (t1597 * t1689 + t1674 * t1745) * t1502;
t1622 = t1343 * t1716;
t1344 = (pkin(6) * t1608 + t1403 * t1483 + pkin(5) + (-t1492 + 0.1e1) * t1742) * t1513 - t1791 * t1671 + (t1596 * t1689 + t1671 * t1745) * t1505;
t1621 = t1344 * t1715;
t1581 = t1682 * t1720;
t1580 = t1682 * t1719;
t1579 = t1682 * t1718;
t1336 = ((-t1695 - m(3) * (rSges(3,2) * t1751 + t1519)) * t1507 + (t1499 * t1444 - m(3) * (rSges(3,1) * t1751 - t1764)) * t1498 - t1418 * t1508) * t1493 - (-Icges(3,5) * t1498 - Icges(3,6) * t1507 + (rSges(3,1) * t1498 + rSges(3,2) * t1507) * m(3) * (rSges(3,3) + t1760)) * t1494;
t1570 = t1336 * t1581;
t1337 = ((-t1694 - m(3) * (rSges(3,2) * t1750 + t1519)) * t1510 + (t1502 * t1444 - m(3) * (rSges(3,1) * t1750 - t1764)) * t1501 - t1418 * t1511) * t1493 - (-Icges(3,5) * t1501 - Icges(3,6) * t1510 + (rSges(3,1) * t1501 + rSges(3,2) * t1510) * m(3) * (rSges(3,3) + t1759)) * t1494;
t1569 = t1337 * t1580;
t1338 = ((-t1693 - m(3) * (rSges(3,2) * t1749 + t1519)) * t1513 + (t1505 * t1444 - m(3) * (rSges(3,1) * t1749 - t1764)) * t1504 - t1418 * t1514) * t1493 - (-Icges(3,5) * t1504 - Icges(3,6) * t1513 + (rSges(3,1) * t1504 + rSges(3,2) * t1513) * m(3) * (rSges(3,3) + t1758)) * t1494;
t1568 = t1338 * t1579;
t1567 = t1375 * t1581;
t1566 = t1418 * t1581;
t1565 = t1376 * t1580;
t1564 = t1418 * t1580;
t1563 = t1377 * t1579;
t1562 = t1418 * t1579;
t1555 = Icges(1,3) + (t1527 + (0.2e1 * pkin(5) + rSges(2,3)) * rSges(2,3) + t1659) * m(2) + (t1525 + t1528 + t1659) * m(3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t1738;
t1541 = -(rSges(2,1) + t1516) * (-rSges(2,1) + t1516) * m(2) + (-t1531 - t1462) * m(3) + t1737 - t1738;
t1467 = 0.2e1 * m(3) * t1519;
t1466 = -0.2e1 * t1662;
t1465 = 0.2e1 * t1662;
t1397 = 0.2e1 * t1405;
t1335 = t1383 * t1470 + t1386 * t1703 + (t1383 * t1665 - t1615) * pkin(2);
t1334 = t1382 * t1469 + t1385 * t1704 + (t1382 * t1667 - t1616) * pkin(2);
t1333 = t1381 * t1468 + t1384 * t1705 + (t1381 * t1669 - t1617) * pkin(2);
t1332 = -t1386 * t1470 + t1383 * t1703 + (-t1386 * t1665 - t1618) * pkin(2);
t1331 = -t1385 * t1469 + t1382 * t1704 + (-t1385 * t1667 - t1619) * pkin(2);
t1330 = -t1384 * t1468 + t1381 * t1705 + (-t1384 * t1669 - t1620) * pkin(2);
t1314 = (0.2e1 * (-t1445 * t1513 + t1696 + t1736) * t1670 + t1397 * t1490 + (t1466 - 0.4e1 * t1672) * t1513 + t1541 + t1809) * t1483 - t1796 * t1654 + ((t1467 + 0.2e1 * t1693) * t1513 + (-t1661 - 0.2e1 * t1696) * t1505 + t1417 * t1771) * t1514 - t1803 + (t1465 + 0.2e1 * t1672) * t1513 + t1505 * t1658 + t1555 - t1809;
t1313 = (0.2e1 * (-t1445 * t1510 + t1697 + t1736) * t1673 + t1397 * t1487 + (t1466 - 0.4e1 * t1675) * t1510 + t1541 + t1810) * t1483 - t1797 * t1654 + ((t1467 + 0.2e1 * t1694) * t1510 + (-t1661 - 0.2e1 * t1697) * t1502 + t1416 * t1771) * t1511 - t1804 + (t1465 + 0.2e1 * t1675) * t1510 + t1502 * t1658 + t1555 - t1810;
t1312 = (0.2e1 * (-t1445 * t1507 + t1698 + t1736) * t1676 + t1397 * t1484 + (t1466 - 0.4e1 * t1678) * t1507 + t1541 + t1811) * t1483 - t1798 * t1654 + ((t1467 + 0.2e1 * t1695) * t1507 + (-t1661 - 0.2e1 * t1698) * t1499 + t1415 * t1771) * t1508 - t1805 + (t1465 + 0.2e1 * t1678) * t1507 + t1499 * t1658 + t1555 - t1811;
t1311 = t1338 * t1712 + (t1341 * t1702 + t1344 * t1377) * t1715;
t1310 = t1337 * t1713 + (t1340 * t1702 + t1343 * t1376) * t1716;
t1309 = t1336 * t1714 + (t1339 * t1702 + t1342 * t1375) * t1717;
t1308 = t1329 * t1712 + (t1344 * t1347 + t1377 * t1724) * t1715;
t1307 = t1328 * t1713 + (t1343 * t1346 + t1376 * t1725) * t1716;
t1306 = t1327 * t1714 + (t1342 * t1345 + t1375 * t1726) * t1717;
t1305 = t1323 * t1562 - (t1320 * t1706 + t1335 * t1338) * t1356;
t1304 = t1322 * t1564 - (t1319 * t1707 + t1334 * t1337) * t1355;
t1303 = t1321 * t1566 - (t1318 * t1708 + t1333 * t1336) * t1354;
t1302 = -t1326 * t1562 - (t1317 * t1706 + t1332 * t1338) * t1356;
t1301 = -t1325 * t1564 - (t1316 * t1707 + t1331 * t1337) * t1355;
t1300 = -t1324 * t1566 - (t1315 * t1708 + t1330 * t1336) * t1354;
t1299 = t1314 * t1712 + (t1329 * t1344 + t1338 * t1724) * t1715;
t1298 = t1313 * t1713 + (t1328 * t1343 + t1337 * t1725) * t1716;
t1297 = t1312 * t1714 + (t1327 * t1342 + t1336 * t1726) * t1717;
t1296 = t1323 * t1563 - (t1320 * t1721 + t1329 * t1335) * t1356;
t1295 = t1322 * t1565 - (t1319 * t1722 + t1328 * t1334) * t1355;
t1294 = t1321 * t1567 - (t1318 * t1723 + t1327 * t1333) * t1354;
t1293 = -t1326 * t1563 - (t1317 * t1721 + t1329 * t1332) * t1356;
t1292 = -t1325 * t1565 - (t1316 * t1722 + t1328 * t1331) * t1355;
t1291 = -t1324 * t1567 - (t1315 * t1723 + t1327 * t1330) * t1354;
t1290 = t1323 * t1568 - (t1314 * t1335 + t1320 * t1727) * t1356;
t1289 = t1322 * t1569 - (t1313 * t1334 + t1319 * t1728) * t1355;
t1288 = t1321 * t1570 - (t1312 * t1333 + t1318 * t1729) * t1354;
t1287 = -t1326 * t1568 - (t1314 * t1332 + t1317 * t1727) * t1356;
t1286 = -t1325 * t1569 - (t1313 * t1331 + t1316 * t1728) * t1355;
t1285 = -t1324 * t1570 - (t1312 * t1330 + t1315 * t1729) * t1354;
t1 = [m(4) - (t1290 * t1335 + t1296 * t1730) * t1356 - (t1289 * t1334 + t1295 * t1731) * t1355 - (t1288 * t1333 + t1294 * t1732) * t1354 + (t1303 * t1632 + t1304 * t1631 + t1305 * t1630) * t1682, -(t1290 * t1332 + t1296 * t1733) * t1356 - (t1289 * t1331 + t1295 * t1734) * t1355 - (t1288 * t1330 + t1294 * t1735) * t1354 + (-t1303 * t1629 - t1304 * t1628 - t1305 * t1627) * t1682, t1294 * t1623 + t1295 * t1622 + t1296 * t1621 + t1288 * t1714 + t1289 * t1713 + t1290 * t1712 + (t1303 * t1626 + t1304 * t1625 + t1305 * t1624) * t1533; -(t1287 * t1335 + t1293 * t1730) * t1356 - (t1286 * t1334 + t1292 * t1731) * t1355 - (t1285 * t1333 + t1291 * t1732) * t1354 + (t1300 * t1632 + t1301 * t1631 + t1302 * t1630) * t1682, m(4) - (t1287 * t1332 + t1293 * t1733) * t1356 - (t1286 * t1331 + t1292 * t1734) * t1355 - (t1285 * t1330 + t1291 * t1735) * t1354 + (-t1300 * t1629 - t1301 * t1628 - t1302 * t1627) * t1682, t1291 * t1623 + t1292 * t1622 + t1293 * t1621 + t1285 * t1714 + t1286 * t1713 + t1287 * t1712 + (t1300 * t1626 + t1301 * t1625 + t1302 * t1624) * t1533; -(t1299 * t1335 + t1308 * t1730) * t1356 - (t1298 * t1334 + t1307 * t1731) * t1355 - (t1297 * t1333 + t1306 * t1732) * t1354 + (t1309 * t1632 + t1310 * t1631 + t1311 * t1630) * t1682, -(t1299 * t1332 + t1308 * t1733) * t1356 - (t1298 * t1331 + t1307 * t1734) * t1355 - (t1297 * t1330 + t1306 * t1735) * t1354 + (-t1309 * t1629 - t1310 * t1628 - t1311 * t1627) * t1682, t1306 * t1623 + t1307 * t1622 + t1308 * t1621 + t1297 * t1714 + t1298 * t1713 + t1299 * t1712 + m(4) + (t1309 * t1626 + t1310 * t1625 + t1311 * t1624) * t1533;];
MX  = t1;