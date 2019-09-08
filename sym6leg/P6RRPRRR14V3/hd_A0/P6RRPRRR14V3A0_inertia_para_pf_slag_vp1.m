% Calculate inertia matrix for parallel robot
% P6RRPRRR14V3G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% MX [6x6]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-15 09:53
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-15 09:47:28
% EndTime: 2019-04-15 09:47:44
% DurationCPUTime: 17.92s
% Computational Cost: add. (33602->1015), mult. (60308->1441), div. (4320->12), fcn. (48836->42), ass. (0->573)
t1769 = xP(4);
t1709 = cos(t1769);
t1767 = xP(6);
t1704 = sin(t1767);
t1707 = cos(t1767);
t1771 = rSges(4,2);
t1772 = rSges(4,1);
t1651 = t1704 * t1772 + t1707 * t1771;
t1706 = sin(t1769);
t2032 = t1651 * t1706;
t1650 = t1704 * t1771 - t1707 * t1772;
t1768 = xP(5);
t1705 = sin(t1768);
t1708 = cos(t1768);
t1770 = rSges(4,3);
t2089 = t1650 * t1705 + t1708 * t1770;
t1560 = m(4) * (t1709 * t2089 + t2032);
t1682 = m(4) * t1770 * t1771 - Icges(4,6);
t2059 = m(4) * t1772;
t1683 = t1770 * t2059 - Icges(4,5);
t1877 = t1704 * t1682 - t1683 * t1707;
t2090 = t1705 * t1877;
t1785 = koppelP(6,3);
t1791 = koppelP(6,2);
t1797 = koppelP(6,1);
t2078 = -t1704 * t1791 + t1707 * t1797;
t1618 = -t1705 * t2078 + t1708 * t1785;
t1658 = t1704 * t1797 + t1707 * t1791;
t1567 = t1618 * t1706 - t1658 * t1709;
t1563 = t1618 * t1709 + t1658 * t1706;
t1786 = koppelP(5,3);
t1792 = koppelP(5,2);
t1798 = koppelP(5,1);
t2079 = -t1704 * t1792 + t1707 * t1798;
t1620 = -t1705 * t2079 + t1708 * t1786;
t1659 = t1704 * t1798 + t1707 * t1792;
t1570 = t1620 * t1706 - t1659 * t1709;
t1564 = t1620 * t1709 + t1659 * t1706;
t1787 = koppelP(4,3);
t1793 = koppelP(4,2);
t1799 = koppelP(4,1);
t2080 = -t1704 * t1793 + t1707 * t1799;
t1622 = -t1705 * t2080 + t1708 * t1787;
t1660 = t1704 * t1799 + t1707 * t1793;
t1573 = t1622 * t1706 - t1660 * t1709;
t1565 = t1622 * t1709 + t1660 * t1706;
t1788 = koppelP(3,3);
t1794 = koppelP(3,2);
t1800 = koppelP(3,1);
t2081 = -t1704 * t1794 + t1707 * t1800;
t1624 = -t1705 * t2081 + t1708 * t1788;
t1661 = t1704 * t1800 + t1707 * t1794;
t1576 = t1624 * t1706 - t1661 * t1709;
t1566 = t1624 * t1709 + t1661 * t1706;
t1789 = koppelP(2,3);
t1795 = koppelP(2,2);
t1801 = koppelP(2,1);
t2082 = -t1704 * t1795 + t1707 * t1801;
t1626 = -t1705 * t2082 + t1708 * t1789;
t1662 = t1704 * t1801 + t1707 * t1795;
t1579 = t1626 * t1706 - t1662 * t1709;
t1580 = t1626 * t1709 + t1662 * t1706;
t1790 = koppelP(1,3);
t1796 = koppelP(1,2);
t1802 = koppelP(1,1);
t2083 = -t1704 * t1796 + t1707 * t1802;
t1628 = -t1705 * t2083 + t1708 * t1790;
t1663 = t1704 * t1802 + t1707 * t1796;
t1583 = t1628 * t1706 - t1663 * t1709;
t1584 = t1628 * t1709 + t1663 * t1706;
t1691 = t1708 ^ 2;
t2087 = -0.2e1 * t1691;
t2086 = m(4) * t1650;
t2060 = m(4) * (t1650 * t1708 - t1705 * t1770);
t1743 = sin(qJ(2,6));
t1749 = cos(qJ(2,6));
t1774 = 0.1e1 / qJ(3,6);
t1737 = legFrame(6,3);
t1692 = sin(t1737);
t1698 = cos(t1737);
t1744 = sin(qJ(1,6));
t1750 = cos(qJ(1,6));
t1639 = t1692 * t1750 + t1698 * t1744;
t2041 = t1563 * t1639;
t1482 = (-t1567 * t1743 - t1749 * t2041) * t1774;
t1488 = t1567 * t1749 - t1743 * t2041;
t2073 = m(3) * rSges(3,2);
t2075 = m(2) * rSges(2,3);
t1671 = rSges(2,1) * t2075 + rSges(3,1) * t2073 - Icges(3,4) - Icges(2,5);
t1888 = -rSges(2,2) * t2075 + Icges(2,6) - Icges(3,6);
t1728 = rSges(3,3) + qJ(3,6);
t2072 = m(3) * t1728;
t1601 = (rSges(3,2) * t2072 + t1888) * t1749 - t1743 * t1671;
t1804 = rSges(3,2) ^ 2;
t1805 = rSges(2,2) ^ 2;
t2046 = -Icges(2,1) - Icges(3,1);
t1886 = Icges(1,3) + (rSges(2,3) ^ 2 + t1805) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t2046;
t1807 = rSges(2,1) ^ 2;
t1887 = (-t1805 + t1807) * m(2) + Icges(2,2) + Icges(3,3) + t2046;
t1889 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t1803 = (rSges(3,3) ^ 2);
t2077 = 2 * rSges(3,3);
t1949 = t1803 + (t2077 + qJ(3,6)) * qJ(3,6);
t1542 = (t1804 + t1949) * m(3) + (0.2e1 * t1743 * (rSges(3,1) * t2072 + t1889) + (-(rSges(3,1) + t1728) * (-rSges(3,1) + t1728) * m(3) + t1887) * t1749) * t1749 + t1886;
t1638 = -t1692 * t1744 + t1698 * t1750;
t2000 = 0.1e1 / t1743 * t1774;
t1962 = t1638 * t2000;
t1919 = t1542 * t1962;
t2058 = rSges(3,2) * t1743;
t1981 = m(3) * t2058;
t1362 = t1482 * t1601 + t1488 * t1981 - t1563 * t1919;
t1745 = sin(qJ(2,5));
t1751 = cos(qJ(2,5));
t1776 = 0.1e1 / qJ(3,5);
t1738 = legFrame(5,3);
t1693 = sin(t1738);
t1699 = cos(t1738);
t1746 = sin(qJ(1,5));
t1752 = cos(qJ(1,5));
t1641 = t1693 * t1752 + t1699 * t1746;
t2040 = t1564 * t1641;
t1483 = (-t1570 * t1745 - t1751 * t2040) * t1776;
t1489 = t1570 * t1751 - t1745 * t2040;
t1729 = rSges(3,3) + qJ(3,5);
t2071 = m(3) * t1729;
t1602 = (rSges(3,2) * t2071 + t1888) * t1751 - t1745 * t1671;
t1948 = t1803 + (t2077 + qJ(3,5)) * qJ(3,5);
t1543 = (t1804 + t1948) * m(3) + (0.2e1 * t1745 * (rSges(3,1) * t2071 + t1889) + (-(rSges(3,1) + t1729) * (-rSges(3,1) + t1729) * m(3) + t1887) * t1751) * t1751 + t1886;
t1640 = -t1693 * t1746 + t1699 * t1752;
t1999 = 0.1e1 / t1745 * t1776;
t1960 = t1640 * t1999;
t1915 = t1543 * t1960;
t2057 = rSges(3,2) * t1745;
t1980 = m(3) * t2057;
t1363 = t1483 * t1602 + t1489 * t1980 - t1564 * t1915;
t1747 = sin(qJ(2,4));
t1753 = cos(qJ(2,4));
t1778 = 0.1e1 / qJ(3,4);
t1739 = legFrame(4,3);
t1694 = sin(t1739);
t1700 = cos(t1739);
t1748 = sin(qJ(1,4));
t1754 = cos(qJ(1,4));
t1643 = t1694 * t1754 + t1700 * t1748;
t2039 = t1565 * t1643;
t1484 = (-t1573 * t1747 - t1753 * t2039) * t1778;
t1490 = t1573 * t1753 - t1747 * t2039;
t1730 = rSges(3,3) + qJ(3,4);
t2070 = m(3) * t1730;
t1603 = (rSges(3,2) * t2070 + t1888) * t1753 - t1747 * t1671;
t1947 = t1803 + (t2077 + qJ(3,4)) * qJ(3,4);
t1544 = (t1804 + t1947) * m(3) + (0.2e1 * t1747 * (rSges(3,1) * t2070 + t1889) + (-(rSges(3,1) + t1730) * (-rSges(3,1) + t1730) * m(3) + t1887) * t1753) * t1753 + t1886;
t1642 = -t1694 * t1748 + t1700 * t1754;
t1998 = 0.1e1 / t1747 * t1778;
t1958 = t1642 * t1998;
t1911 = t1544 * t1958;
t2056 = rSges(3,2) * t1747;
t1979 = m(3) * t2056;
t1364 = t1484 * t1603 + t1490 * t1979 - t1565 * t1911;
t1755 = sin(qJ(2,3));
t1761 = cos(qJ(2,3));
t1780 = 0.1e1 / qJ(3,3);
t1740 = legFrame(3,3);
t1695 = sin(t1740);
t1701 = cos(t1740);
t1756 = sin(qJ(1,3));
t1762 = cos(qJ(1,3));
t1645 = t1695 * t1762 + t1701 * t1756;
t2038 = t1566 * t1645;
t1485 = (-t1576 * t1755 - t1761 * t2038) * t1780;
t1491 = t1576 * t1761 - t1755 * t2038;
t1731 = rSges(3,3) + qJ(3,3);
t2069 = m(3) * t1731;
t1604 = (rSges(3,2) * t2069 + t1888) * t1761 - t1755 * t1671;
t1946 = t1803 + (t2077 + qJ(3,3)) * qJ(3,3);
t1545 = (t1804 + t1946) * m(3) + (0.2e1 * t1755 * (rSges(3,1) * t2069 + t1889) + (-(rSges(3,1) + t1731) * (-rSges(3,1) + t1731) * m(3) + t1887) * t1761) * t1761 + t1886;
t1644 = -t1695 * t1756 + t1701 * t1762;
t1997 = 0.1e1 / t1755 * t1780;
t1956 = t1644 * t1997;
t1907 = t1545 * t1956;
t2055 = rSges(3,2) * t1755;
t1978 = m(3) * t2055;
t1365 = t1485 * t1604 + t1491 * t1978 - t1566 * t1907;
t1757 = sin(qJ(2,2));
t1763 = cos(qJ(2,2));
t1782 = 0.1e1 / qJ(3,2);
t1741 = legFrame(2,3);
t1696 = sin(t1741);
t1702 = cos(t1741);
t1758 = sin(qJ(1,2));
t1764 = cos(qJ(1,2));
t1647 = t1696 * t1764 + t1702 * t1758;
t2036 = t1580 * t1647;
t1486 = (-t1579 * t1757 - t1763 * t2036) * t1782;
t1492 = t1579 * t1763 - t1757 * t2036;
t1732 = rSges(3,3) + qJ(3,2);
t2068 = m(3) * t1732;
t1605 = (rSges(3,2) * t2068 + t1888) * t1763 - t1757 * t1671;
t1945 = t1803 + (t2077 + qJ(3,2)) * qJ(3,2);
t1546 = (t1804 + t1945) * m(3) + (0.2e1 * t1757 * (rSges(3,1) * t2068 + t1889) + (-(rSges(3,1) + t1732) * (-rSges(3,1) + t1732) * m(3) + t1887) * t1763) * t1763 + t1886;
t1646 = -t1696 * t1758 + t1702 * t1764;
t1996 = 0.1e1 / t1757 * t1782;
t1954 = t1646 * t1996;
t1903 = t1546 * t1954;
t2054 = rSges(3,2) * t1757;
t1977 = m(3) * t2054;
t1366 = t1486 * t1605 + t1492 * t1977 - t1580 * t1903;
t1759 = sin(qJ(2,1));
t1765 = cos(qJ(2,1));
t1784 = 0.1e1 / qJ(3,1);
t1742 = legFrame(1,3);
t1697 = sin(t1742);
t1703 = cos(t1742);
t1760 = sin(qJ(1,1));
t1766 = cos(qJ(1,1));
t1649 = t1697 * t1766 + t1703 * t1760;
t2034 = t1584 * t1649;
t1487 = (-t1583 * t1759 - t1765 * t2034) * t1784;
t1493 = t1583 * t1765 - t1759 * t2034;
t1733 = rSges(3,3) + qJ(3,1);
t2067 = m(3) * t1733;
t1606 = (rSges(3,2) * t2067 + t1888) * t1765 - t1759 * t1671;
t1944 = t1803 + (t2077 + qJ(3,1)) * qJ(3,1);
t1547 = (t1804 + t1944) * m(3) + (0.2e1 * t1759 * (rSges(3,1) * t2067 + t1889) + (-(rSges(3,1) + t1733) * (-rSges(3,1) + t1733) * m(3) + t1887) * t1765) * t1765 + t1886;
t1648 = -t1697 * t1760 + t1703 * t1766;
t1995 = 0.1e1 / t1759 * t1784;
t1952 = t1648 * t1995;
t1899 = t1547 * t1952;
t2053 = rSges(3,2) * t1759;
t1976 = m(3) * t2053;
t1367 = t1487 * t1606 + t1493 * t1976 - t1584 * t1899;
t1806 = rSges(3,1) ^ 2;
t1963 = Icges(3,2) + Icges(2,3) + (t1805 + t1807) * m(2);
t1632 = (t1806 + t1949) * m(3) + t1963;
t1918 = t1601 * t1962;
t2074 = m(3) * rSges(3,1);
t1386 = t1482 * t1632 - t1488 * t2074 - t1563 * t1918;
t1633 = (t1806 + t1948) * m(3) + t1963;
t1914 = t1602 * t1960;
t1387 = t1483 * t1633 - t1489 * t2074 - t1564 * t1914;
t1634 = (t1806 + t1947) * m(3) + t1963;
t1910 = t1603 * t1958;
t1388 = t1484 * t1634 - t1490 * t2074 - t1565 * t1910;
t1635 = (t1806 + t1946) * m(3) + t1963;
t1906 = t1604 * t1956;
t1389 = t1485 * t1635 - t1491 * t2074 - t1566 * t1906;
t1636 = (t1806 + t1945) * m(3) + t1963;
t1902 = t1605 * t1954;
t1390 = t1486 * t1636 - t1492 * t2074 - t1580 * t1902;
t1637 = (t1806 + t1944) * m(3) + t1963;
t1898 = t1606 * t1952;
t1391 = t1487 * t1637 - t1493 * t2074 - t1584 * t1898;
t2052 = rSges(3,2) * t1774;
t1975 = t1638 * t2052;
t1410 = (-rSges(3,1) * t1482 - t1563 * t1975 + t1488) * m(3);
t2051 = rSges(3,2) * t1776;
t1973 = t1640 * t2051;
t1411 = (-rSges(3,1) * t1483 - t1564 * t1973 + t1489) * m(3);
t2050 = rSges(3,2) * t1778;
t1971 = t1642 * t2050;
t1412 = (-rSges(3,1) * t1484 - t1565 * t1971 + t1490) * m(3);
t2049 = rSges(3,2) * t1780;
t1969 = t1644 * t2049;
t1413 = (-rSges(3,1) * t1485 - t1566 * t1969 + t1491) * m(3);
t2048 = rSges(3,2) * t1782;
t1967 = t1646 * t2048;
t1414 = (-rSges(3,1) * t1486 - t1580 * t1967 + t1492) * m(3);
t2047 = rSges(3,2) * t1784;
t1965 = t1648 * t2047;
t1415 = (-rSges(3,1) * t1487 - t1584 * t1965 + t1493) * m(3);
t1735 = t1771 ^ 2;
t1736 = t1772 ^ 2;
t1679 = m(4) * (t1735 + t1736) + Icges(4,3);
t1932 = t1584 * t1952;
t1934 = t1580 * t1954;
t1936 = t1566 * t1956;
t1937 = t1565 * t1958;
t1938 = t1564 * t1960;
t1939 = t1563 * t1962;
t2084 = t1679 - t1362 * t1939 - t1363 * t1938 - t1364 * t1937 - t1365 * t1936 - t1366 * t1934 - t1367 * t1932 + t1386 * t1482 + t1387 * t1483 + t1388 * t1484 + t1389 * t1485 + t1390 * t1486 + t1391 * t1487 + t1410 * t1488 + t1411 * t1489 + t1412 * t1490 + t1413 * t1491 + t1414 * t1492 + t1415 * t1493;
t1684 = t1771 * t2059 - Icges(4,4);
t2076 = -0.2e1 * t1684;
t2066 = m(3) * t1749;
t2065 = m(3) * t1751;
t2064 = m(3) * t1753;
t2063 = m(3) * t1761;
t2062 = m(3) * t1763;
t2061 = m(3) * t1765;
t2045 = t1563 * t1638;
t2044 = t1564 * t1640;
t2043 = t1565 * t1642;
t2042 = t1566 * t1644;
t2037 = t1580 * t1646;
t2035 = t1584 * t1648;
t2031 = t1651 * t1709;
t2018 = t1683 * t1704;
t2017 = t1684 * t1705;
t2016 = t1704 * t1707;
t2009 = t1705 * t1706;
t2008 = t1705 * t1709;
t2001 = t1709 * t1684;
t1994 = t1743 * t1774;
t1993 = t1745 * t1776;
t1992 = t1747 * t1778;
t1991 = t1749 * t1774;
t1990 = t1751 * t1776;
t1989 = t1753 * t1778;
t1988 = t1755 * t1780;
t1987 = t1757 * t1782;
t1986 = t1759 * t1784;
t1985 = t1761 * t1780;
t1984 = t1763 * t1782;
t1983 = t1765 * t1784;
t1734 = t1770 ^ 2;
t1982 = (t1734 - t1735) * m(4) + Icges(4,2) - Icges(4,3);
t1974 = t1639 * t2052;
t1972 = t1641 * t2051;
t1970 = t1643 * t2050;
t1968 = t1645 * t2049;
t1966 = t1647 * t2048;
t1964 = t1649 * t2047;
t1961 = t1639 * t2000;
t1959 = t1641 * t1999;
t1957 = t1643 * t1998;
t1955 = t1645 * t1997;
t1953 = t1647 * t1996;
t1951 = t1649 * t1995;
t1950 = t1684 * t2016;
t1943 = t1563 * t1961;
t1942 = t1564 * t1959;
t1941 = t1565 * t1957;
t1940 = t1566 * t1955;
t1935 = t1580 * t1953;
t1933 = t1584 * t1951;
t1588 = (-rSges(3,2) * t2066 + t1601 * t1774) * t1743;
t1931 = t1588 * t1962;
t1930 = t1588 * t1961;
t1589 = (-rSges(3,2) * t2065 + t1602 * t1776) * t1745;
t1929 = t1589 * t1960;
t1928 = t1589 * t1959;
t1590 = (-rSges(3,2) * t2064 + t1603 * t1778) * t1747;
t1927 = t1590 * t1958;
t1926 = t1590 * t1957;
t1591 = (-rSges(3,2) * t2063 + t1604 * t1780) * t1755;
t1925 = t1591 * t1956;
t1924 = t1591 * t1955;
t1592 = (-rSges(3,2) * t2062 + t1605 * t1782) * t1757;
t1923 = t1592 * t1954;
t1922 = t1592 * t1953;
t1593 = (-rSges(3,2) * t2061 + t1606 * t1784) * t1759;
t1921 = t1593 * t1952;
t1920 = t1593 * t1951;
t1917 = t1542 * t1961;
t1916 = t1601 * t1961;
t1913 = t1543 * t1959;
t1912 = t1602 * t1959;
t1909 = t1544 * t1957;
t1908 = t1603 * t1957;
t1905 = t1545 * t1955;
t1904 = t1604 * t1955;
t1901 = t1546 * t1953;
t1900 = t1605 * t1953;
t1897 = t1547 * t1951;
t1896 = t1606 * t1951;
t1895 = -rSges(3,1) * t1991 + t1743;
t1894 = -rSges(3,1) * t1990 + t1745;
t1893 = -rSges(3,1) * t1989 + t1747;
t1892 = -rSges(3,1) * t1985 + t1755;
t1891 = -rSges(3,1) * t1984 + t1757;
t1890 = -rSges(3,1) * t1983 + t1759;
t1612 = t1705 * t1785 + t1708 * t2078;
t1885 = t1567 * t1638 + t1612 * t1639;
t1613 = t1705 * t1786 + t1708 * t2079;
t1884 = t1570 * t1640 + t1613 * t1641;
t1614 = t1705 * t1787 + t1708 * t2080;
t1883 = t1573 * t1642 + t1614 * t1643;
t1615 = t1705 * t1788 + t1708 * t2081;
t1882 = t1576 * t1644 + t1615 * t1645;
t1616 = t1705 * t1789 + t1708 * t2082;
t1881 = t1579 * t1646 + t1616 * t1647;
t1617 = t1705 * t1790 + t1708 * t2083;
t1880 = t1583 * t1648 + t1617 * t1649;
t1878 = t1682 * t1707 + t2018;
t1494 = t1885 * t1991;
t1495 = (-t1567 * t1639 + t1612 * t1638) * t2000;
t1506 = t1885 * t1743;
t1332 = t1494 * t1601 + t1495 * t1542 + t1506 * t1981;
t1496 = t1884 * t1990;
t1497 = (-t1570 * t1641 + t1613 * t1640) * t1999;
t1507 = t1884 * t1745;
t1333 = t1496 * t1602 + t1497 * t1543 + t1507 * t1980;
t1498 = t1883 * t1989;
t1499 = (-t1573 * t1643 + t1614 * t1642) * t1998;
t1508 = t1883 * t1747;
t1334 = t1498 * t1603 + t1499 * t1544 + t1508 * t1979;
t1500 = t1882 * t1985;
t1501 = (-t1576 * t1645 + t1615 * t1644) * t1997;
t1509 = t1882 * t1755;
t1335 = t1500 * t1604 + t1501 * t1545 + t1509 * t1978;
t1502 = t1881 * t1984;
t1503 = (-t1579 * t1647 + t1616 * t1646) * t1996;
t1510 = t1881 * t1757;
t1336 = t1502 * t1605 + t1503 * t1546 + t1510 * t1977;
t1504 = t1880 * t1983;
t1505 = (-t1583 * t1649 + t1617 * t1648) * t1995;
t1511 = t1880 * t1759;
t1337 = t1504 * t1606 + t1505 * t1547 + t1511 * t1976;
t1368 = t1494 * t1632 + t1495 * t1601 - t1506 * t2074;
t1369 = t1496 * t1633 + t1497 * t1602 - t1507 * t2074;
t1370 = t1498 * t1634 + t1499 * t1603 - t1508 * t2074;
t1371 = t1500 * t1635 + t1501 * t1604 - t1509 * t2074;
t1372 = t1502 * t1636 + t1503 * t1605 - t1510 * t2074;
t1373 = t1504 * t1637 + t1505 * t1606 - t1511 * t2074;
t1392 = (-rSges(3,1) * t1494 + t1495 * t2058 + t1506) * m(3);
t1393 = (-rSges(3,1) * t1496 + t1497 * t2057 + t1507) * m(3);
t1394 = (-rSges(3,1) * t1498 + t1499 * t2056 + t1508) * m(3);
t1398 = (-rSges(3,1) * t1500 + t1501 * t2055 + t1509) * m(3);
t1399 = (-rSges(3,1) * t1502 + t1503 * t2054 + t1510) * m(3);
t1400 = (-rSges(3,1) * t1504 + t1505 * t2053 + t1511) * m(3);
t1876 = t1332 * t1495 + t1333 * t1497 + t1334 * t1499 + t1335 * t1501 + t1336 * t1503 + t1368 * t1494 + t1369 * t1496 + t1370 * t1498 + t1371 * t1500 + t1372 * t1502 + t1392 * t1506 + t1393 * t1507 + t1394 * t1508 + t1398 * t1509 + t1399 * t1510 + t1337 * t1505 + t1373 * t1504 + t1400 * t1511;
t1530 = (-t1612 * t1743 + t1749 * t2045) * t1774;
t1531 = (-t1613 * t1745 + t1751 * t2044) * t1776;
t1532 = (-t1614 * t1747 + t1753 * t2043) * t1778;
t1533 = (-t1615 * t1755 + t1761 * t2042) * t1780;
t1534 = (-t1616 * t1757 + t1763 * t2037) * t1782;
t1535 = (-t1617 * t1759 + t1765 * t2035) * t1784;
t1536 = t1612 * t1749 + t1743 * t2045;
t1537 = t1613 * t1751 + t1745 * t2044;
t1538 = t1614 * t1753 + t1747 * t2043;
t1539 = t1615 * t1761 + t1755 * t2042;
t1540 = t1616 * t1763 + t1757 * t2037;
t1541 = t1617 * t1765 + t1759 * t2035;
t1875 = t1332 * t1943 + t1333 * t1942 + t1334 * t1941 + t1335 * t1940 + t1336 * t1935 - t1368 * t1530 - t1369 * t1531 - t1370 * t1532 - t1371 * t1533 - t1372 * t1534 - t1392 * t1536 - t1393 * t1537 - t1394 * t1538 - t1398 * t1539 - t1399 * t1540 + t1337 * t1933 - t1373 * t1535 - t1400 * t1541;
t1874 = t1362 * t1495 + t1363 * t1497 + t1364 * t1499 + t1365 * t1501 + t1366 * t1503 + t1386 * t1494 + t1387 * t1496 + t1388 * t1498 + t1389 * t1500 + t1390 * t1502 + t1410 * t1506 + t1411 * t1507 + t1412 * t1508 + t1413 * t1509 + t1414 * t1510 + t1367 * t1505 + t1391 * t1504 + t1415 * t1511;
t1873 = -t1362 * t1943 - t1363 * t1942 - t1364 * t1941 - t1365 * t1940 - t1366 * t1935 + t1386 * t1530 + t1387 * t1531 + t1388 * t1532 + t1389 * t1533 + t1390 * t1534 + t1410 * t1536 + t1411 * t1537 + t1412 * t1538 + t1413 * t1539 + t1414 * t1540 - t1367 * t1933 + t1391 * t1535 + t1415 * t1541;
t1380 = t1530 * t1601 + t1536 * t1981 - t1563 * t1917;
t1381 = t1531 * t1602 + t1537 * t1980 - t1564 * t1913;
t1382 = t1532 * t1603 + t1538 * t1979 - t1565 * t1909;
t1383 = t1533 * t1604 + t1539 * t1978 - t1566 * t1905;
t1384 = t1534 * t1605 + t1540 * t1977 - t1580 * t1901;
t1385 = t1535 * t1606 + t1541 * t1976 - t1584 * t1897;
t1404 = t1530 * t1632 - t1536 * t2074 - t1563 * t1916;
t1405 = t1531 * t1633 - t1537 * t2074 - t1564 * t1912;
t1406 = t1532 * t1634 - t1538 * t2074 - t1565 * t1908;
t1407 = t1533 * t1635 - t1539 * t2074 - t1566 * t1904;
t1408 = t1534 * t1636 - t1540 * t2074 - t1580 * t1900;
t1409 = t1535 * t1637 - t1541 * t2074 - t1584 * t1896;
t1416 = (-rSges(3,1) * t1530 - t1563 * t1974 + t1536) * m(3);
t1417 = (-rSges(3,1) * t1531 - t1564 * t1972 + t1537) * m(3);
t1418 = (-rSges(3,1) * t1532 - t1565 * t1970 + t1538) * m(3);
t1419 = (-rSges(3,1) * t1533 - t1566 * t1968 + t1539) * m(3);
t1420 = (-rSges(3,1) * t1534 - t1580 * t1966 + t1540) * m(3);
t1421 = (-rSges(3,1) * t1535 - t1584 * t1964 + t1541) * m(3);
t1872 = t1380 * t1495 + t1381 * t1497 + t1382 * t1499 + t1383 * t1501 + t1384 * t1503 + t1404 * t1494 + t1405 * t1496 + t1406 * t1498 + t1407 * t1500 + t1408 * t1502 + t1416 * t1506 + t1417 * t1507 + t1418 * t1508 + t1419 * t1509 + t1420 * t1510 + t1385 * t1505 + t1409 * t1504 + t1421 * t1511;
t1871 = t1380 * t1943 + t1381 * t1942 + t1382 * t1941 + t1383 * t1940 + t1384 * t1935 - t1404 * t1530 - t1405 * t1531 - t1406 * t1532 - t1407 * t1533 - t1408 * t1534 - t1416 * t1536 - t1417 * t1537 - t1418 * t1538 - t1419 * t1539 - t1420 * t1540 + t1385 * t1933 - t1409 * t1535 - t1421 * t1541;
t1820 = t1743 ^ 2 * t2073 + t1601 * t1991;
t1458 = t1639 * t1820 + t1919;
t1819 = t1745 ^ 2 * t2073 + t1602 * t1990;
t1460 = t1641 * t1819 + t1915;
t1818 = t1747 ^ 2 * t2073 + t1603 * t1989;
t1462 = t1643 * t1818 + t1911;
t1817 = t1755 ^ 2 * t2073 + t1604 * t1985;
t1464 = t1645 * t1817 + t1907;
t1816 = t1757 ^ 2 * t2073 + t1605 * t1984;
t1466 = t1647 * t1816 + t1903;
t1815 = t1759 ^ 2 * t2073 + t1606 * t1983;
t1468 = t1649 * t1815 + t1899;
t1826 = t1632 * t1991 - t1743 * t2074;
t1512 = t1639 * t1826 + t1918;
t1825 = t1633 * t1990 - t1745 * t2074;
t1514 = t1641 * t1825 + t1914;
t1824 = t1634 * t1989 - t1747 * t2074;
t1516 = t1643 * t1824 + t1910;
t1823 = t1635 * t1985 - t1755 * t2074;
t1521 = t1645 * t1823 + t1906;
t1822 = t1636 * t1984 - t1757 * t2074;
t1523 = t1647 * t1822 + t1902;
t1821 = t1637 * t1983 - t1759 * t2074;
t1525 = t1649 * t1821 + t1898;
t1548 = (t1639 * t1895 + t1975) * m(3);
t1550 = (t1641 * t1894 + t1973) * m(3);
t1552 = (t1643 * t1893 + t1971) * m(3);
t1554 = (t1645 * t1892 + t1969) * m(3);
t1556 = (t1647 * t1891 + t1967) * m(3);
t1558 = (t1649 * t1890 + t1965) * m(3);
t1870 = t1458 * t1495 + t1460 * t1497 + t1462 * t1499 + t1464 * t1501 + t1466 * t1503 + t1494 * t1512 + t1496 * t1514 + t1498 * t1516 + t1500 * t1521 + t1502 * t1523 + t1506 * t1548 + t1507 * t1550 + t1508 * t1552 + t1509 * t1554 + t1510 * t1556 + t1468 * t1505 + t1504 * t1525 + t1511 * t1558;
t1459 = t1638 * t1820 - t1917;
t1461 = t1640 * t1819 - t1913;
t1463 = t1642 * t1818 - t1909;
t1465 = t1644 * t1817 - t1905;
t1467 = t1646 * t1816 - t1901;
t1469 = t1648 * t1815 - t1897;
t1513 = t1638 * t1826 - t1916;
t1515 = t1640 * t1825 - t1912;
t1517 = t1642 * t1824 - t1908;
t1522 = t1644 * t1823 - t1904;
t1524 = t1646 * t1822 - t1900;
t1526 = t1648 * t1821 - t1896;
t1549 = (t1638 * t1895 - t1974) * m(3);
t1551 = (t1640 * t1894 - t1972) * m(3);
t1553 = (t1642 * t1893 - t1970) * m(3);
t1555 = (t1644 * t1892 - t1968) * m(3);
t1557 = (t1646 * t1891 - t1966) * m(3);
t1559 = (t1648 * t1890 - t1964) * m(3);
t1869 = t1459 * t1495 + t1461 * t1497 + t1463 * t1499 + t1465 * t1501 + t1467 * t1503 + t1494 * t1513 + t1496 * t1515 + t1498 * t1517 + t1500 * t1522 + t1502 * t1524 + t1506 * t1549 + t1507 * t1551 + t1508 * t1553 + t1509 * t1555 + t1510 * t1557 + t1469 * t1505 + t1504 * t1526 + t1511 * t1559;
t1868 = -t1458 * t1943 - t1460 * t1942 - t1462 * t1941 - t1464 * t1940 - t1466 * t1935 + t1512 * t1530 + t1514 * t1531 + t1516 * t1532 + t1521 * t1533 + t1523 * t1534 + t1536 * t1548 + t1537 * t1550 + t1538 * t1552 + t1539 * t1554 + t1540 * t1556 - t1468 * t1933 + t1525 * t1535 + t1541 * t1558;
t1867 = -t1459 * t1943 - t1461 * t1942 - t1463 * t1941 - t1465 * t1940 - t1467 * t1935 + t1513 * t1530 + t1515 * t1531 + t1517 * t1532 + t1522 * t1533 + t1524 * t1534 + t1536 * t1549 + t1537 * t1551 + t1538 * t1553 + t1539 * t1555 + t1540 * t1557 - t1469 * t1933 + t1526 * t1535 + t1541 * t1559;
t1598 = rSges(3,1) * t2066 + t1632 * t1994;
t1599 = rSges(3,1) * t2065 + t1633 * t1993;
t1600 = rSges(3,1) * t2064 + t1634 * t1992;
t1607 = rSges(3,1) * t2063 + t1635 * t1988;
t1608 = rSges(3,1) * t2062 + t1636 * t1987;
t1609 = rSges(3,1) * t2061 + t1637 * t1986;
t1664 = (-rSges(3,1) * t1994 - t1749) * m(3);
t1665 = (-rSges(3,1) * t1993 - t1751) * m(3);
t1666 = (-rSges(3,1) * t1992 - t1753) * m(3);
t1667 = (-rSges(3,1) * t1988 - t1761) * m(3);
t1668 = (-rSges(3,1) * t1987 - t1763) * m(3);
t1669 = (-rSges(3,1) * t1986 - t1765) * m(3);
t1866 = t1494 * t1598 + t1495 * t1588 + t1496 * t1599 + t1497 * t1589 + t1498 * t1600 + t1499 * t1590 + t1500 * t1607 + t1501 * t1591 + t1502 * t1608 + t1503 * t1592 + t1506 * t1664 + t1507 * t1665 + t1508 * t1666 + t1509 * t1667 + t1510 * t1668 + t1504 * t1609 + t1505 * t1593 + t1511 * t1669;
t1865 = -t1530 * t1598 - t1531 * t1599 - t1532 * t1600 - t1533 * t1607 - t1534 * t1608 - t1536 * t1664 - t1537 * t1665 - t1538 * t1666 - t1539 * t1667 - t1540 * t1668 + t1563 * t1930 + t1564 * t1928 + t1565 * t1926 + t1566 * t1924 + t1580 * t1922 - t1535 * t1609 - t1541 * t1669 + t1584 * t1920;
t1864 = t1368 * t1991 + t1392 * t1743;
t1863 = t1369 * t1990 + t1393 * t1745;
t1862 = t1370 * t1989 + t1394 * t1747;
t1861 = t1371 * t1985 + t1398 * t1755;
t1860 = t1372 * t1984 + t1399 * t1757;
t1859 = t1373 * t1983 + t1400 * t1759;
t1858 = t1386 * t1991 + t1410 * t1743;
t1857 = t1387 * t1990 + t1411 * t1745;
t1856 = t1388 * t1989 + t1412 * t1747;
t1855 = t1389 * t1985 + t1413 * t1755;
t1854 = t1390 * t1984 + t1414 * t1757;
t1853 = t1391 * t1983 + t1415 * t1759;
t1852 = t1404 * t1991 + t1416 * t1743;
t1851 = t1405 * t1990 + t1417 * t1745;
t1850 = t1406 * t1989 + t1418 * t1747;
t1849 = t1407 * t1985 + t1419 * t1755;
t1848 = t1408 * t1984 + t1420 * t1757;
t1847 = t1409 * t1983 + t1421 * t1759;
t1846 = t1512 * t1991 + t1548 * t1743;
t1845 = t1513 * t1991 + t1549 * t1743;
t1844 = t1514 * t1990 + t1550 * t1745;
t1843 = t1515 * t1990 + t1551 * t1745;
t1842 = t1516 * t1989 + t1552 * t1747;
t1841 = t1517 * t1989 + t1553 * t1747;
t1840 = t1521 * t1985 + t1554 * t1755;
t1839 = t1522 * t1985 + t1555 * t1755;
t1838 = t1523 * t1984 + t1556 * t1757;
t1837 = t1524 * t1984 + t1557 * t1757;
t1836 = t1525 * t1983 + t1558 * t1759;
t1835 = t1526 * t1983 + t1559 * t1759;
t1670 = (t1735 - t1736) * m(4) + Icges(4,1) - Icges(4,2);
t1834 = -(t1670 * t1704 * t1705 + t1682 * t1708) * t1707 - t1708 * t2018;
t1833 = t1598 * t1991 + t1664 * t1743;
t1832 = t1599 * t1990 + t1665 * t1745;
t1831 = t1600 * t1989 + t1666 * t1747;
t1830 = t1607 * t1985 + t1667 * t1755;
t1829 = t1608 * t1984 + t1668 * t1757;
t1828 = t1609 * t1983 + t1669 * t1759;
t1690 = t1707 ^ 2;
t1827 = t1670 * t1690 + 0.2e1 * t1950;
t1813 = -t1332 * t1939 - t1333 * t1938 - t1334 * t1937 - t1335 * t1936 - t1336 * t1934 - t1337 * t1932 + t1368 * t1482 + t1369 * t1483 + t1370 * t1484 + t1371 * t1485 + t1372 * t1486 + t1373 * t1487 + t1392 * t1488 + t1393 * t1489 + t1394 * t1490 + t1398 * t1491 + t1399 * t1492 + t1400 * t1493;
t1811 = -t1380 * t1939 - t1381 * t1938 - t1382 * t1937 - t1383 * t1936 - t1384 * t1934 - t1385 * t1932 + t1404 * t1482 + t1405 * t1483 + t1406 * t1484 + t1407 * t1485 + t1408 * t1486 + t1409 * t1487 + t1416 * t1488 + t1417 * t1489 + t1418 * t1490 + t1419 * t1491 + t1420 * t1492 + t1421 * t1493;
t1810 = -t1458 * t1939 - t1460 * t1938 - t1462 * t1937 - t1464 * t1936 - t1466 * t1934 - t1468 * t1932 + t1512 * t1482 + t1514 * t1483 + t1516 * t1484 + t1521 * t1485 + t1523 * t1486 + t1525 * t1487 + t1548 * t1488 + t1550 * t1489 + t1552 * t1490 + t1554 * t1491 + t1556 * t1492 + t1558 * t1493;
t1809 = -t1459 * t1939 - t1461 * t1938 - t1463 * t1937 - t1465 * t1936 - t1467 * t1934 - t1469 * t1932 + t1513 * t1482 + t1515 * t1483 + t1517 * t1484 + t1522 * t1485 + t1524 * t1486 + t1526 * t1487 + t1549 * t1488 + t1551 * t1489 + t1553 * t1490 + t1555 * t1491 + t1557 * t1492 + t1559 * t1493;
t1808 = t1598 * t1482 + t1599 * t1483 + t1600 * t1484 + t1607 * t1485 + t1608 * t1486 + t1609 * t1487 + t1664 * t1488 + t1665 * t1489 + t1666 * t1490 + t1667 * t1491 + t1668 * t1492 + t1669 * t1493 - t1563 * t1931 - t1564 * t1929 - t1565 * t1927 - t1566 * t1925 - t1580 * t1923 - t1584 * t1921;
t1678 = m(4) * (t1734 + t1735) + Icges(4,1);
t1596 = -t1670 * t1709 + t2009 * t2076;
t1595 = t1670 * t2016 + t1690 * t2076 + t1684;
t1594 = t1827 + t1982;
t1587 = -t1708 * t1679 + t2090;
t1561 = m(4) * (-t1706 * t2089 + t2031);
t1 = [-t1459 * t1961 - t1461 * t1959 - t1463 * t1957 - t1465 * t1955 - t1467 * t1953 - t1469 * t1951 + t1638 * t1845 + t1640 * t1843 + t1642 * t1841 + t1644 * t1839 + t1646 * t1837 + t1648 * t1835 + m(4), t1459 * t1962 + t1461 * t1960 + t1463 * t1958 + t1465 * t1956 + t1467 * t1954 + t1469 * t1952 + t1639 * t1845 + t1641 * t1843 + t1643 * t1841 + t1645 * t1839 + t1647 * t1837 + t1649 * t1835, t1513 * t1994 + t1515 * t1993 + t1517 * t1992 + t1522 * t1988 + t1524 * t1987 + t1526 * t1986 - t1549 * t1749 - t1551 * t1751 - t1553 * t1753 - t1555 * t1761 - t1557 * t1763 - t1559 * t1765, t1809, m(4) * t2089 + t1869 * t1706 + t1867 * t1709, t1809 * t1705 + (-m(4) * t1651 - t1706 * t1867 + t1709 * t1869) * t1708; -t1458 * t1961 - t1460 * t1959 - t1462 * t1957 - t1464 * t1955 - t1466 * t1953 - t1468 * t1951 + t1638 * t1846 + t1640 * t1844 + t1642 * t1842 + t1644 * t1840 + t1646 * t1838 + t1648 * t1836, t1458 * t1962 + t1460 * t1960 + t1462 * t1958 + t1464 * t1956 + t1466 * t1954 + t1468 * t1952 + t1639 * t1846 + t1641 * t1844 + t1643 * t1842 + t1645 * t1840 + t1647 * t1838 + t1649 * t1836 + m(4), t1512 * t1994 + t1514 * t1993 + t1516 * t1992 + t1521 * t1988 + t1523 * t1987 + t1525 * t1986 - t1548 * t1749 - t1550 * t1751 - t1552 * t1753 - t1554 * t1761 - t1556 * t1763 - t1558 * t1765, t1810 - t1560, t1868 * t1709 + (t1870 - t2060) * t1706, -t1709 * t2086 + (-m(4) * t2032 + t1810) * t1705 + (-t1706 * t1868 + t1709 * t1870) * t1708; t1638 * t1833 + t1640 * t1832 + t1642 * t1831 + t1644 * t1830 + t1646 * t1829 + t1648 * t1828 - t1920 - t1922 - t1924 - t1926 - t1928 - t1930, t1639 * t1833 + t1641 * t1832 + t1643 * t1831 + t1645 * t1830 + t1647 * t1829 + t1649 * t1828 + t1921 + t1923 + t1925 + t1927 + t1929 + t1931, t1598 * t1994 + t1599 * t1993 + t1600 * t1992 + t1607 * t1988 + t1608 * t1987 + t1609 * t1986 - t1664 * t1749 - t1665 * t1751 - t1666 * t1753 - t1667 * t1761 - t1668 * t1763 - t1669 * t1765 + m(4), t1808 + t1561, t1866 * t1706 + (-t1865 + t2060) * t1709, -t1706 * t2086 + (m(4) * t2031 + t1808) * t1705 + (t1706 * t1865 + t1709 * t1866) * t1708; -t1362 * t1961 - t1363 * t1959 - t1364 * t1957 - t1365 * t1955 - t1366 * t1953 - t1367 * t1951 + t1638 * t1858 + t1640 * t1857 + t1642 * t1856 + t1644 * t1855 + t1646 * t1854 + t1648 * t1853, t1362 * t1962 + t1363 * t1960 + t1364 * t1958 + t1365 * t1956 + t1366 * t1954 + t1367 * t1952 + t1639 * t1858 + t1641 * t1857 + t1643 * t1856 + t1645 * t1855 + t1647 * t1854 + t1649 * t1853 - t1560, t1386 * t1994 + t1387 * t1993 + t1388 * t1992 + t1389 * t1988 + t1390 * t1987 + t1391 * t1986 - t1410 * t1749 - t1411 * t1751 - t1412 * t1753 - t1413 * t1761 - t1414 * t1763 - t1415 * t1765 + t1561, t1594 * t1691 + 0.2e1 * t1708 * t2090 + t2084, t1595 * t1708 - t1705 * t1878 + t1706 * t1874 + t1709 * t1873, t2084 * t1705 + (-t1706 * t1873 + t1709 * t1874 + t1877) * t1708; -t1380 * t1961 - t1381 * t1959 - t1382 * t1957 - t1383 * t1955 - t1384 * t1953 - t1385 * t1951 + t1852 * t1638 + t1851 * t1640 + t1850 * t1642 + t1849 * t1644 + t1848 * t1646 + t1847 * t1648 + t1560, t1380 * t1962 + t1381 * t1960 + t1382 * t1958 + t1383 * t1956 + t1384 * t1954 + t1385 * t1952 + t1639 * t1852 + t1641 * t1851 + t1643 * t1850 + t1645 * t1849 + t1647 * t1848 + t1649 * t1847, t1404 * t1994 + t1405 * t1993 + t1406 * t1992 + t1407 * t1988 + t1408 * t1987 + t1409 * t1986 - t1416 * t1749 - t1417 * t1751 - t1418 * t1753 - t1419 * t1761 - t1420 * t1763 - t1421 * t1765 + t2060 ((t1670 * t2009 - 0.2e1 * t2001) * t1690 - t1596 * t2016 + t1982 * t2009 + t2001) * t1708 + (-t1682 * t2008 - t1683 * t1706) * t1707 - t1704 * (-t1682 * t1706 + t1683 * t2008) + t1811 + t1706 * t1877 * t2087, t1596 * t1690 + (t1678 - t1871 - 0.2e1 * t1950) * t1709 + (-t1834 + t1872 + t2017) * t1706, t1587 * t1706 - t1878 * t1709 + t1811 * t1705 + (t1706 * t1871 + t1709 * t1872) * t1708; -t1332 * t1961 - t1333 * t1959 - t1334 * t1957 - t1335 * t1955 - t1336 * t1953 - t1337 * t1951 + t1638 * t1864 + t1640 * t1863 + t1642 * t1862 + t1644 * t1861 + t1646 * t1860 + t1648 * t1859 - t1561, t1332 * t1962 + t1333 * t1960 + t1334 * t1958 + t1335 * t1956 + t1336 * t1954 + t1337 * t1952 + t1639 * t1864 + t1641 * t1863 + t1643 * t1862 + t1645 * t1861 + t1647 * t1860 + t1649 * t1859 - t2060, t1368 * t1994 + t1369 * t1993 + t1370 * t1992 + t1371 * t1988 + t1372 * t1987 + t1373 * t1986 - t1392 * t1749 - t1393 * t1751 - t1394 * t1753 - t1398 * t1761 - t1399 * t1763 - t1400 * t1765 (-t1594 * t2008 + t1595 * t1706) * t1708 - (t2087 + 0.1e1) * t1709 * t1877 + t1813 - t1878 * t2009 ((0.2e1 * t1690 - 0.1e1) * t2017 + t1834 - t1875) * t1709 + (t1678 - t1827 + t1876) * t1706, -t1587 * t1709 - t1706 * t1878 + t1813 * t1705 + (t1706 * t1875 + t1709 * t1876) * t1708;];
MX  = t1;
