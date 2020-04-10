% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G2P3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G2P3A0_coriolisvec_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:43
% EndTime: 2020-03-09 21:16:46
% DurationCPUTime: 3.14s
% Computational Cost: add. (1746->237), mult. (5668->560), div. (2193->17), fcn. (5130->18), ass. (0->256)
t1682 = cos(qJ(3,3));
t1655 = 0.1e1 / t1682;
t1704 = t1682 ^ 2;
t1656 = 0.1e1 / t1704;
t1657 = t1655 * t1656;
t1690 = 0.1e1 / pkin(2);
t1673 = legFrame(3,2);
t1634 = sin(t1673);
t1637 = cos(t1673);
t1689 = xDP(1);
t1879 = xDP(2);
t1622 = -t1634 * t1879 + t1637 * t1689;
t1619 = t1622 ^ 2;
t1677 = sin(qJ(2,3));
t1642 = 0.1e1 / t1677;
t1860 = t1619 * t1642;
t1688 = xDP(3);
t1676 = sin(qJ(3,3));
t1683 = cos(qJ(2,3));
t1826 = t1676 * t1683;
t1610 = t1622 * t1826 - t1682 * t1688;
t1607 = t1610 ^ 2;
t1643 = 0.1e1 / t1677 ^ 2;
t1869 = t1607 * t1643;
t1868 = t1642 * t1869;
t1598 = (t1860 + t1868) * t1690 * t1657;
t1641 = t1676 ^ 2;
t1660 = t1683 ^ 2;
t1847 = t1642 * t1660;
t1777 = t1656 * t1847;
t1896 = t1598 * (t1641 * t1777 - t1677);
t1895 = t1598 * (t1677 + t1847);
t1684 = cos(qJ(3,2));
t1661 = 0.1e1 / t1684;
t1708 = t1684 ^ 2;
t1662 = 0.1e1 / t1708;
t1663 = t1661 * t1662;
t1674 = legFrame(2,2);
t1635 = sin(t1674);
t1638 = cos(t1674);
t1624 = -t1635 * t1879 + t1638 * t1689;
t1620 = t1624 ^ 2;
t1679 = sin(qJ(2,2));
t1647 = 0.1e1 / t1679;
t1857 = t1620 * t1647;
t1678 = sin(qJ(3,2));
t1685 = cos(qJ(2,2));
t1823 = t1678 * t1685;
t1611 = t1624 * t1823 - t1684 * t1688;
t1608 = t1611 ^ 2;
t1648 = 0.1e1 / t1679 ^ 2;
t1867 = t1608 * t1648;
t1866 = t1647 * t1867;
t1599 = (t1857 + t1866) * t1690 * t1663;
t1646 = t1678 ^ 2;
t1666 = t1685 ^ 2;
t1843 = t1647 * t1666;
t1775 = t1662 * t1843;
t1894 = t1599 * (t1646 * t1775 - t1679);
t1893 = t1599 * (t1679 + t1843);
t1686 = cos(qJ(3,1));
t1667 = 0.1e1 / t1686;
t1712 = t1686 ^ 2;
t1668 = 0.1e1 / t1712;
t1669 = t1667 * t1668;
t1675 = legFrame(1,2);
t1636 = sin(t1675);
t1639 = cos(t1675);
t1626 = -t1636 * t1879 + t1639 * t1689;
t1621 = t1626 ^ 2;
t1681 = sin(qJ(2,1));
t1652 = 0.1e1 / t1681;
t1854 = t1621 * t1652;
t1680 = sin(qJ(3,1));
t1687 = cos(qJ(2,1));
t1820 = t1680 * t1687;
t1612 = t1626 * t1820 - t1686 * t1688;
t1609 = t1612 ^ 2;
t1653 = 0.1e1 / t1681 ^ 2;
t1865 = t1609 * t1653;
t1864 = t1652 * t1865;
t1600 = (t1854 + t1864) * t1690 * t1669;
t1651 = t1680 ^ 2;
t1672 = t1687 ^ 2;
t1839 = t1652 * t1672;
t1773 = t1668 * t1839;
t1892 = t1600 * (t1651 * t1773 - t1681);
t1891 = t1600 * (t1681 + t1839);
t1658 = 0.1e1 / t1704 ^ 2;
t1691 = 0.1e1 / pkin(2) ^ 2;
t1801 = t1691 * t1869;
t1604 = t1658 * t1801;
t1827 = t1676 * t1677;
t1846 = t1642 * t1683;
t1851 = t1622 * t1683;
t1863 = t1610 * t1643;
t1880 = t1690 ^ 2;
t1586 = (-(t1610 * t1846 + t1622 * t1827) * t1863 + (-t1610 * t1676 - t1851) * t1642 * t1622) * t1657 * t1655 * t1880;
t1810 = t1586 * t1846;
t1890 = t1655 * (0.2e1 * t1641 * t1810 + 0.2e1 * t1656 * t1801 - t1604);
t1664 = 0.1e1 / t1708 ^ 2;
t1799 = t1691 * t1867;
t1605 = t1664 * t1799;
t1824 = t1678 * t1679;
t1842 = t1647 * t1685;
t1850 = t1624 * t1685;
t1862 = t1611 * t1648;
t1587 = (-(t1611 * t1842 + t1624 * t1824) * t1862 + (-t1611 * t1678 - t1850) * t1647 * t1624) * t1663 * t1661 * t1880;
t1808 = t1587 * t1842;
t1889 = t1661 * (0.2e1 * t1646 * t1808 + 0.2e1 * t1662 * t1799 - t1605);
t1670 = 0.1e1 / t1712 ^ 2;
t1797 = t1691 * t1865;
t1606 = t1670 * t1797;
t1821 = t1680 * t1681;
t1838 = t1652 * t1687;
t1849 = t1626 * t1687;
t1861 = t1612 * t1653;
t1588 = (-(t1612 * t1838 + t1626 * t1821) * t1861 + (-t1612 * t1680 - t1849) * t1652 * t1626) * t1669 * t1667 * t1880;
t1806 = t1588 * t1838;
t1888 = t1667 * (0.2e1 * t1651 * t1806 + 0.2e1 * t1668 * t1797 - t1606);
t1770 = 0.2e1 * t1612 * t1849;
t1837 = t1653 * t1670;
t1887 = (-t1609 * t1680 + t1651 * t1770) * t1837;
t1771 = 0.2e1 * t1611 * t1850;
t1841 = t1648 * t1664;
t1886 = (-t1608 * t1678 + t1646 * t1771) * t1841;
t1772 = 0.2e1 * t1610 * t1851;
t1845 = t1643 * t1658;
t1885 = (-t1607 * t1676 + t1641 * t1772) * t1845;
t1650 = t1680 * t1651;
t1671 = t1667 * t1670;
t1829 = t1669 * t1680;
t1884 = (t1650 * t1671 + t1829) * t1621 * t1838;
t1645 = t1678 * t1646;
t1665 = t1661 * t1664;
t1832 = t1663 * t1678;
t1883 = (t1645 * t1665 + t1832) * t1620 * t1842;
t1640 = t1676 * t1641;
t1659 = t1655 * t1658;
t1835 = t1657 * t1676;
t1882 = (t1640 * t1659 + t1835) * t1619 * t1846;
t1881 = 2 * MDP(6);
t1878 = MDP(2) * t1690;
t1877 = MDP(8) * t1690;
t1692 = t1690 * t1691;
t1876 = MDP(9) * t1692;
t1875 = t1586 * t1655;
t1874 = t1586 * t1683;
t1873 = t1587 * t1661;
t1872 = t1587 * t1685;
t1871 = t1588 * t1667;
t1870 = t1588 * t1687;
t1859 = t1619 * t1658;
t1858 = t1619 * t1691;
t1856 = t1620 * t1664;
t1855 = t1620 * t1691;
t1853 = t1621 * t1670;
t1852 = t1621 * t1691;
t1848 = t1642 * t1655;
t1844 = t1647 * t1661;
t1840 = t1652 * t1667;
t1836 = t1655 * t1676;
t1834 = t1659 * t1683;
t1833 = t1661 * t1678;
t1831 = t1665 * t1685;
t1830 = t1667 * t1680;
t1828 = t1671 * t1687;
t1825 = t1677 * t1682;
t1822 = t1679 * t1684;
t1819 = t1681 * t1686;
t1601 = t1656 * t1858 + t1604;
t1730 = t1642 * t1691 * t1772;
t1791 = t1677 * t1858;
t1818 = -t1641 * t1657 * t1791 - t1601 * t1825 + t1682 * t1874 + t1730 * t1835;
t1602 = t1662 * t1855 + t1605;
t1729 = t1647 * t1691 * t1771;
t1788 = t1679 * t1855;
t1817 = -t1646 * t1663 * t1788 - t1602 * t1822 + t1684 * t1872 + t1729 * t1832;
t1603 = t1668 * t1852 + t1606;
t1728 = t1652 * t1691 * t1770;
t1785 = t1681 * t1852;
t1816 = -t1651 * t1669 * t1785 - t1603 * t1819 + t1686 * t1870 + t1728 * t1829;
t1815 = (-t1656 * t1791 - t1874) * t1676 + t1601 * t1827 + t1656 * t1730;
t1814 = (-t1662 * t1788 - t1872) * t1678 + t1602 * t1824 + t1662 * t1729;
t1813 = (-t1668 * t1785 - t1870) * t1680 + t1603 * t1821 + t1668 * t1728;
t1812 = 0.2e1 * t1692;
t1811 = t1586 * t1848;
t1809 = t1587 * t1844;
t1807 = t1588 * t1840;
t1805 = t1598 * t1848;
t1804 = t1599 * t1844;
t1803 = t1600 * t1840;
t1802 = t1659 * t1869;
t1800 = t1665 * t1867;
t1798 = t1671 * t1865;
t1796 = t1622 * t1863;
t1795 = t1624 * t1862;
t1794 = t1626 * t1861;
t1792 = t1676 * t1859;
t1789 = t1678 * t1856;
t1786 = t1680 * t1853;
t1784 = t1634 * t1836;
t1783 = t1635 * t1833;
t1782 = t1636 * t1830;
t1781 = t1637 * t1836;
t1780 = t1638 * t1833;
t1779 = t1639 * t1830;
t1778 = t1655 * t1846;
t1776 = t1661 * t1842;
t1774 = t1667 * t1838;
t1766 = t1586 * t1778;
t1765 = t1656 * t1810;
t1764 = t1587 * t1776;
t1763 = t1662 * t1808;
t1762 = t1588 * t1774;
t1761 = t1668 * t1806;
t1760 = t1598 * t1778;
t1759 = t1598 * t1656 * t1826;
t1758 = t1599 * t1776;
t1757 = t1599 * t1662 * t1823;
t1756 = t1600 * t1774;
t1755 = t1600 * t1668 * t1820;
t1754 = t1834 * t1868;
t1753 = t1831 * t1866;
t1752 = t1828 * t1864;
t1751 = t1676 * t1796;
t1750 = t1678 * t1795;
t1749 = t1680 * t1794;
t1748 = t1818 * t1848;
t1747 = t1817 * t1844;
t1746 = t1816 * t1840;
t1745 = t1815 * t1848;
t1744 = t1814 * t1844;
t1743 = t1813 * t1840;
t1742 = t1640 * t1765;
t1741 = t1676 * t1765;
t1740 = t1645 * t1763;
t1739 = t1678 * t1763;
t1738 = t1650 * t1761;
t1737 = t1680 * t1761;
t1736 = t1598 * t1676 * t1777;
t1735 = t1599 * t1678 * t1775;
t1734 = t1600 * t1680 * t1773;
t1631 = -0.1e1 + 0.2e1 * t1704;
t1721 = t1631 * t1751 * t1834;
t1632 = -0.1e1 + 0.2e1 * t1708;
t1720 = t1632 * t1750 * t1831;
t1633 = -0.1e1 + 0.2e1 * t1712;
t1719 = t1633 * t1749 * t1828;
t1618 = t1636 * t1680 + t1639 * t1819;
t1617 = t1635 * t1678 + t1638 * t1822;
t1616 = t1634 * t1676 + t1637 * t1825;
t1615 = t1636 * t1819 - t1639 * t1680;
t1614 = t1635 * t1822 - t1638 * t1678;
t1613 = t1634 * t1825 - t1637 * t1676;
t1 = [(t1613 * t1805 + t1614 * t1804 + t1615 * t1803) * MDP(1) + (t1637 * t1741 + t1638 * t1739 + t1639 * t1737) * t1878 + (t1613 * t1766 + t1614 * t1764 + t1615 * t1762 + (-t1613 * t1802 - t1614 * t1800 - t1615 * t1798) * t1691 + (t1637 * t1736 + t1638 * t1735 + t1639 * t1734) * t1690) * MDP(3) + (-t1613 * t1875 - t1614 * t1873 - t1615 * t1871 + (-t1613 * t1754 - t1614 * t1753 - t1615 * t1752) * t1691 + (-t1637 * t1759 - t1638 * t1757 - t1639 * t1755) * t1690) * MDP(4) + ((t1637 * t1742 + t1638 * t1740 + t1639 * t1738) * t1690 + (-t1637 * t1885 - t1638 * t1886 - t1639 * t1887) * t1692) * MDP(5) + ((-t1637 * t1721 - t1638 * t1720 - t1639 * t1719) * t1812 + (t1637 * t1890 + t1638 * t1889 + t1639 * t1888) * t1690) * MDP(6) + ((-t1586 * t1781 - t1587 * t1780 - t1588 * t1779) * t1690 + (t1637 * t1882 + t1638 * t1883 + t1639 * t1884) * t1692) * MDP(7) + (-t1586 * t1637 - t1587 * t1638 - t1588 * t1639) * t1877 + (-t1637 * t1792 - t1638 * t1789 - t1639 * t1786) * t1876 + (t1615 * t1746 + t1614 * t1747 + t1613 * t1748 + (t1779 * t1891 + t1780 * t1893 + t1781 * t1895) * t1690) * MDP(10) + (t1615 * t1743 + t1614 * t1744 + t1613 * t1745 + (-t1637 * t1896 - t1638 * t1894 - t1639 * t1892) * t1690) * MDP(11); (t1616 * t1805 + t1617 * t1804 + t1618 * t1803) * MDP(1) + (-t1634 * t1741 - t1635 * t1739 - t1636 * t1737) * t1878 + (t1616 * t1766 + t1617 * t1764 + t1618 * t1762 + (-t1616 * t1802 - t1617 * t1800 - t1618 * t1798) * t1691 + (-t1634 * t1736 - t1635 * t1735 - t1636 * t1734) * t1690) * MDP(3) + (-t1616 * t1875 - t1617 * t1873 - t1618 * t1871 + (-t1616 * t1754 - t1617 * t1753 - t1618 * t1752) * t1691 + (t1634 * t1759 + t1635 * t1757 + t1636 * t1755) * t1690) * MDP(4) + ((-t1634 * t1742 - t1635 * t1740 - t1636 * t1738) * t1690 + (t1634 * t1885 + t1635 * t1886 + t1636 * t1887) * t1692) * MDP(5) + ((t1634 * t1721 + t1635 * t1720 + t1636 * t1719) * t1812 + (-t1634 * t1890 - t1635 * t1889 - t1636 * t1888) * t1690) * MDP(6) + ((t1586 * t1784 + t1587 * t1783 + t1588 * t1782) * t1690 + (-t1634 * t1882 - t1635 * t1883 - t1636 * t1884) * t1692) * MDP(7) + (t1586 * t1634 + t1587 * t1635 + t1588 * t1636) * t1877 + (t1634 * t1792 + t1635 * t1789 + t1636 * t1786) * t1876 + (t1618 * t1746 + t1617 * t1747 + t1616 * t1748 + (-t1782 * t1891 - t1783 * t1893 - t1784 * t1895) * t1690) * MDP(10) + (t1618 * t1743 + t1617 * t1744 + t1616 * t1745 + (t1634 * t1896 + t1635 * t1894 + t1636 * t1892) * t1690) * MDP(11); (-t1870 - t1872 - t1874) * MDP(4) + (t1588 * t1672 * MDP(3) + (t1600 * MDP(1) + MDP(10) * t1816 + MDP(11) * t1813) * t1687) * t1652 + (t1587 * t1666 * MDP(3) + (t1599 * MDP(1) + MDP(10) * t1817 + MDP(11) * t1814) * t1685) * t1647 + (t1586 * t1660 * MDP(3) + (t1598 * MDP(1) + MDP(10) * t1818 + MDP(11) * t1815) * t1683) * t1642 + ((-t1607 * t1683 * t1845 - t1608 * t1685 * t1841 - t1609 * t1687 * t1837) * MDP(3) + (-t1658 * t1660 * t1868 - t1664 * t1666 * t1866 - t1670 * t1672 * t1864) * MDP(4)) * t1691 + (0.2e1 * (t1657 * t1751 + t1663 * t1750 + t1669 * t1749) * MDP(5) + (t1631 * t1658 * t1796 + t1632 * t1664 * t1795 + t1633 * t1670 * t1794) * t1881 + (-t1641 * t1642 * t1859 - t1646 * t1647 * t1856 - t1651 * t1652 * t1853 - t1656 * t1860 - t1662 * t1857 - t1668 * t1854) * MDP(7)) * t1692 + ((-t1807 - t1809 - t1811) * MDP(2) + (-t1756 - t1758 - t1760) * MDP(3) + (t1598 * t1655 + t1599 * t1661 + t1600 * t1667) * MDP(4) + (-t1641 * t1811 - t1646 * t1809 - t1651 * t1807) * MDP(5) + (-t1586 * t1642 * t1676 - t1587 * t1647 * t1678 - t1588 * t1652 * t1680) * t1881 + (-t1598 * t1846 - t1599 * t1842 - t1600 * t1838) * MDP(10) + (t1676 * t1760 + t1678 * t1758 + t1680 * t1756) * MDP(11)) * t1690;];
taucX  = t1;
