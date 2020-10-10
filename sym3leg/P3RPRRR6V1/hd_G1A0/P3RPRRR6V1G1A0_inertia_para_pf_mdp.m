% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR6V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRRR6V1G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:32
% EndTime: 2020-08-06 18:32:34
% DurationCPUTime: 2.41s
% Computational Cost: add. (7135->284), mult. (5491->489), div. (516->14), fcn. (3096->65), ass. (0->241)
t1934 = 2 * pkin(2);
t1933 = MDP(10) / 0.2e1;
t1932 = MDP(11) / 0.2e1;
t1798 = pkin(1) ^ 2;
t1931 = MDP(4) * t1798 + MDP(1);
t1784 = sin(qJ(3,3));
t1787 = cos(qJ(3,3));
t1930 = t1784 * t1787;
t1785 = sin(qJ(3,2));
t1788 = cos(qJ(3,2));
t1929 = t1785 * t1788;
t1786 = sin(qJ(3,1));
t1789 = cos(qJ(3,1));
t1928 = t1786 * t1789;
t1927 = -2 * pkin(1);
t1925 = 2 * MDP(6);
t1774 = qJ(1,3) + legFrame(3,3);
t1763 = pkin(7) + t1774;
t1745 = sin(t1763);
t1924 = -0.2e1 * t1745;
t1923 = 0.2e1 * t1745;
t1775 = qJ(1,2) + legFrame(2,3);
t1764 = pkin(7) + t1775;
t1746 = sin(t1764);
t1922 = -0.2e1 * t1746;
t1921 = 0.2e1 * t1746;
t1776 = qJ(1,1) + legFrame(1,3);
t1765 = pkin(7) + t1776;
t1747 = sin(t1765);
t1920 = -0.2e1 * t1747;
t1919 = 0.2e1 * t1747;
t1748 = cos(t1763);
t1918 = -0.2e1 * t1748;
t1917 = 0.2e1 * t1748;
t1749 = cos(t1764);
t1916 = -0.2e1 * t1749;
t1915 = 0.2e1 * t1749;
t1750 = cos(t1765);
t1914 = -0.2e1 * t1750;
t1913 = 0.2e1 * t1750;
t1792 = (pkin(6) + pkin(5));
t1912 = -2 * t1792;
t1911 = 2 * t1792;
t1773 = cos(pkin(7)) * pkin(1) + pkin(2);
t1910 = MDP(4) / 0.2e1;
t1909 = MDP(4) / 0.4e1;
t1907 = MDP(9) / pkin(3) ^ 2;
t1906 = 2 * pkin(1);
t1752 = qJ(3,3) + t1763;
t1753 = -qJ(3,3) + t1763;
t1864 = sin(t1752) + sin(t1753);
t1703 = t1748 * t1911 + sin(t1774) * t1927 + pkin(2) * t1924 - t1864 * pkin(3);
t1793 = 0.2e1 * qJ(3,3);
t1715 = pkin(3) * sin(t1793) + t1784 * t1934 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1);
t1709 = 0.1e1 / t1715;
t1904 = t1703 * t1709;
t1756 = qJ(3,2) + t1764;
t1757 = -qJ(3,2) + t1764;
t1863 = sin(t1756) + sin(t1757);
t1704 = t1749 * t1911 + sin(t1775) * t1927 + pkin(2) * t1922 - t1863 * pkin(3);
t1794 = 0.2e1 * qJ(3,2);
t1716 = pkin(3) * sin(t1794) + t1785 * t1934 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1);
t1711 = 0.1e1 / t1716;
t1903 = t1704 * t1711;
t1760 = qJ(3,1) + t1765;
t1761 = -qJ(3,1) + t1765;
t1862 = sin(t1760) + sin(t1761);
t1705 = t1750 * t1911 + sin(t1776) * t1927 + pkin(2) * t1920 - t1862 * pkin(3);
t1795 = 0.2e1 * qJ(3,1);
t1717 = pkin(3) * sin(t1795) + t1786 * t1934 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1);
t1713 = 0.1e1 / t1717;
t1902 = t1705 * t1713;
t1861 = cos(t1752) + cos(t1753);
t1706 = t1745 * t1912 + cos(t1774) * t1927 + pkin(2) * t1918 - t1861 * pkin(3);
t1901 = t1706 * t1709;
t1860 = cos(t1756) + cos(t1757);
t1707 = t1746 * t1912 + cos(t1775) * t1927 + pkin(2) * t1916 - t1860 * pkin(3);
t1900 = t1707 * t1711;
t1859 = cos(t1760) + cos(t1761);
t1708 = t1747 * t1912 + cos(t1776) * t1927 + pkin(2) * t1914 - t1859 * pkin(3);
t1899 = t1708 * t1713;
t1898 = t1709 * t1784;
t1897 = t1709 * t1787;
t1896 = t1709 / 0.2e1;
t1710 = 0.1e1 / t1715 ^ 2;
t1895 = t1710 * t1784;
t1894 = t1710 * t1787;
t1893 = t1711 * t1785;
t1892 = t1711 * t1788;
t1891 = t1711 / 0.2e1;
t1712 = 0.1e1 / t1716 ^ 2;
t1890 = t1712 * t1785;
t1889 = t1712 * t1788;
t1888 = t1713 * t1786;
t1887 = t1713 * t1789;
t1886 = t1713 / 0.2e1;
t1714 = 0.1e1 / t1717 ^ 2;
t1885 = t1714 * t1786;
t1884 = t1714 * t1789;
t1724 = t1787 * pkin(3) + t1773;
t1718 = 0.1e1 / t1724;
t1883 = t1718 * t1745;
t1882 = t1718 * t1748;
t1881 = t1718 * t1784;
t1880 = t1718 * t1787;
t1725 = t1788 * pkin(3) + t1773;
t1720 = 0.1e1 / t1725;
t1879 = t1720 * t1746;
t1878 = t1720 * t1749;
t1877 = t1720 * t1785;
t1876 = t1720 * t1788;
t1726 = t1789 * pkin(3) + t1773;
t1722 = 0.1e1 / t1726;
t1875 = t1722 * t1747;
t1874 = t1722 * t1750;
t1873 = t1722 * t1786;
t1872 = t1722 * t1789;
t1719 = 0.1e1 / t1724 ^ 2;
t1739 = t1745 ^ 2;
t1871 = t1739 * t1719;
t1721 = 0.1e1 / t1725 ^ 2;
t1740 = t1746 ^ 2;
t1870 = t1740 * t1721;
t1723 = 0.1e1 / t1726 ^ 2;
t1741 = t1747 ^ 2;
t1869 = t1741 * t1723;
t1742 = t1748 ^ 2;
t1868 = t1742 * t1719;
t1743 = t1749 ^ 2;
t1867 = t1743 * t1721;
t1744 = t1750 ^ 2;
t1866 = t1744 * t1723;
t1772 = sin(pkin(7)) * pkin(1) + pkin(5);
t1796 = 0.1e1 / pkin(3);
t1865 = t1772 * t1796;
t1858 = 2 * MDP(7);
t1857 = 2 * MDP(8);
t1856 = t1703 * t1895;
t1855 = t1703 * t1894;
t1854 = t1704 * t1890;
t1853 = t1704 * t1889;
t1852 = t1705 * t1885;
t1851 = t1705 * t1884;
t1850 = t1706 * t1895;
t1849 = t1706 * t1894;
t1848 = t1707 * t1890;
t1847 = t1707 * t1889;
t1846 = t1708 * t1885;
t1845 = t1708 * t1884;
t1844 = t1709 * t1865;
t1843 = t1784 * t1896;
t1842 = t1787 * t1896;
t1841 = t1711 * t1865;
t1840 = t1785 * t1891;
t1839 = t1788 * t1891;
t1838 = t1713 * t1865;
t1837 = t1786 * t1886;
t1836 = t1789 * t1886;
t1835 = t1772 * t1881;
t1834 = t1772 * t1880;
t1833 = t1773 * t1881;
t1832 = t1773 * t1880;
t1831 = t1719 * t1930;
t1830 = t1772 * t1877;
t1829 = t1772 * t1876;
t1828 = t1773 * t1877;
t1827 = t1773 * t1876;
t1826 = t1721 * t1929;
t1825 = t1772 * t1873;
t1824 = t1772 * t1872;
t1823 = t1773 * t1873;
t1822 = t1773 * t1872;
t1821 = t1723 * t1928;
t1820 = t1745 * t1719 * t1748;
t1819 = t1746 * t1721 * t1749;
t1818 = t1747 * t1723 * t1750;
t1817 = t1882 * t1904;
t1816 = t1878 * t1903;
t1815 = t1874 * t1902;
t1814 = t1883 * t1901;
t1813 = t1879 * t1900;
t1812 = t1875 * t1899;
t1811 = t1784 * t1844;
t1810 = t1787 * t1844;
t1809 = t1785 * t1841;
t1808 = t1788 * t1841;
t1807 = t1786 * t1838;
t1806 = t1789 * t1838;
t1751 = t1793 + t1763;
t1754 = -0.2e1 * qJ(3,3) + t1763;
t1766 = qJ(3,3) + t1774;
t1767 = -qJ(3,3) + t1774;
t1695 = t1861 * t1934 + (cos(t1767) + cos(t1766)) * t1906 + t1864 * t1911 + (cos(t1754) + cos(t1751) + t1917) * pkin(3);
t1755 = t1794 + t1764;
t1758 = -0.2e1 * qJ(3,2) + t1764;
t1768 = qJ(3,2) + t1775;
t1769 = -qJ(3,2) + t1775;
t1696 = t1860 * t1934 + (cos(t1769) + cos(t1768)) * t1906 + t1863 * t1911 + (cos(t1758) + cos(t1755) + t1915) * pkin(3);
t1759 = t1795 + t1765;
t1762 = -0.2e1 * qJ(3,1) + t1765;
t1770 = qJ(3,1) + t1776;
t1771 = -qJ(3,1) + t1776;
t1697 = t1859 * t1934 + (cos(t1771) + cos(t1770)) * t1906 + t1862 * t1911 + (cos(t1762) + cos(t1759) + t1913) * pkin(3);
t1698 = t1864 * t1934 + (sin(t1767) + sin(t1766)) * t1906 + t1861 * t1912 + (sin(t1754) + sin(t1751) + t1923) * pkin(3);
t1699 = t1863 * t1934 + (sin(t1769) + sin(t1768)) * t1906 + t1860 * t1912 + (sin(t1758) + sin(t1755) + t1921) * pkin(3);
t1700 = t1862 * t1934 + (sin(t1771) + sin(t1770)) * t1906 + t1859 * t1912 + (sin(t1762) + sin(t1759) + t1919) * pkin(3);
t1778 = t1784 ^ 2;
t1779 = t1785 ^ 2;
t1780 = t1786 ^ 2;
t1799 = -t1818 - t1819 - t1820;
t1802 = t1713 * t1722 * (-t1705 * t1747 + t1708 * t1750);
t1803 = t1711 * t1720 * (-t1704 * t1746 + t1707 * t1749);
t1804 = t1709 * t1718 * (-t1703 * t1745 + t1706 * t1748);
t1805 = (t1799 * t1798 + t1695 * t1698 * t1710 / 0.4e1 + t1696 * t1699 * t1712 / 0.4e1 + t1697 * t1700 * t1714 / 0.4e1) * MDP(4) + (t1703 * t1706 * t1710 + t1704 * t1707 * t1712 + t1705 * t1708 * t1714) * t1907 + (-t1818 * t1928 - t1819 * t1929 - t1820 * t1930) * t1925 + (-t1778 * t1820 - t1779 * t1819 - t1780 * t1818) * MDP(5) + t1799 * MDP(1) + ((t1784 * t1804 + t1785 * t1803 + t1786 * t1802) * MDP(7) + (t1787 * t1804 + t1788 * t1803 + t1789 * t1802) * MDP(8)) * t1796;
t1693 = -t1705 * t1807 + t1822 * t1913;
t1692 = -t1708 * t1807 + t1822 * t1920;
t1691 = -t1703 * t1811 + t1832 * t1917;
t1690 = -t1706 * t1811 + t1832 * t1924;
t1689 = -t1704 * t1809 + t1827 * t1915;
t1688 = -t1707 * t1809 + t1827 * t1922;
t1687 = -t1708 * t1806 + t1823 * t1919;
t1686 = -t1705 * t1806 + t1823 * t1914;
t1685 = -t1707 * t1808 + t1828 * t1921;
t1684 = -t1704 * t1808 + t1828 * t1916;
t1683 = -t1706 * t1810 + t1833 * t1923;
t1682 = -t1703 * t1810 + t1833 * t1918;
t1681 = -t1700 * t1837 - t1750 * t1824;
t1680 = -t1697 * t1837 + t1747 * t1824;
t1679 = -t1699 * t1840 - t1749 * t1829;
t1678 = t1699 * t1839 - t1749 * t1830;
t1677 = -t1698 * t1843 - t1748 * t1834;
t1676 = t1698 * t1842 - t1748 * t1835;
t1675 = t1700 * t1836 - t1750 * t1825;
t1674 = t1697 * t1836 + t1747 * t1825;
t1673 = -t1696 * t1840 + t1746 * t1829;
t1672 = t1696 * t1839 + t1746 * t1830;
t1671 = -t1695 * t1843 + t1745 * t1834;
t1670 = t1695 * t1842 + t1745 * t1835;
t1665 = (t1698 * t1709 + t1699 * t1711 + t1700 * t1713) * t1910 + ((t1703 * t1897 + t1704 * t1892 + t1705 * t1887) * MDP(10) + (-t1703 * t1898 - t1704 * t1893 - t1705 * t1888) * MDP(11)) * t1796;
t1664 = (t1695 * t1709 + t1696 * t1711 + t1697 * t1713) * t1910 + ((t1706 * t1897 + t1707 * t1892 + t1708 * t1887) * MDP(10) + (-t1706 * t1898 - t1707 * t1893 - t1708 * t1888) * MDP(11)) * t1796;
t1 = [(t1778 * t1871 + t1779 * t1870 + t1780 * t1869) * MDP(5) + (t1739 * t1831 + t1740 * t1826 + t1741 * t1821) * t1925 + (-t1688 * t1879 - t1690 * t1883 - t1692 * t1875) * MDP(10) + (-t1683 * t1883 - t1685 * t1879 - t1687 * t1875) * MDP(11) + MDP(12) + (t1695 ^ 2 * t1710 + t1696 ^ 2 * t1712 + t1697 ^ 2 * t1714) * t1909 + (t1706 ^ 2 * t1710 + t1707 ^ 2 * t1712 + t1708 ^ 2 * t1714) * t1907 + ((-t1784 * t1814 - t1785 * t1813 - t1786 * t1812) * t1858 + (-t1787 * t1814 - t1788 * t1813 - t1789 * t1812) * t1857 + (t1670 * t1901 + t1672 * t1900 + t1674 * t1899) * MDP(10) + (t1671 * t1901 + t1673 * t1900 + t1680 * t1899) * MDP(11) + (t1695 * t1849 + t1696 * t1847 + t1697 * t1845) * t1933 + (-t1695 * t1850 - t1696 * t1848 - t1697 * t1846) * t1932) * t1796 + t1931 * (t1869 + t1870 + t1871); (-t1689 * t1879 - t1691 * t1883 - t1693 * t1875) * MDP(10) + (-t1682 * t1883 - t1684 * t1879 - t1686 * t1875) * MDP(11) + ((t1675 * t1899 + t1676 * t1901 + t1678 * t1900) * MDP(10) + (t1677 * t1901 + t1679 * t1900 + t1681 * t1899) * MDP(11) + (t1695 * t1855 + t1696 * t1853 + t1697 * t1851) * t1933 + (-t1695 * t1856 - t1696 * t1854 - t1697 * t1852) * t1932) * t1796 + t1805; t1664; (t1688 * t1878 + t1690 * t1882 + t1692 * t1874) * MDP(10) + (t1683 * t1882 + t1685 * t1878 + t1687 * t1874) * MDP(11) + ((t1670 * t1904 + t1672 * t1903 + t1674 * t1902) * MDP(10) + (t1671 * t1904 + t1673 * t1903 + t1680 * t1902) * MDP(11) + (t1698 * t1849 + t1699 * t1847 + t1700 * t1845) * t1933 + (-t1698 * t1850 - t1699 * t1848 - t1700 * t1846) * t1932) * t1796 + t1805; (t1778 * t1868 + t1779 * t1867 + t1780 * t1866) * MDP(5) + (t1742 * t1831 + t1743 * t1826 + t1744 * t1821) * t1925 + (t1689 * t1878 + t1691 * t1882 + t1693 * t1874) * MDP(10) + (t1682 * t1882 + t1684 * t1878 + t1686 * t1874) * MDP(11) + MDP(12) + (t1698 ^ 2 * t1710 + t1699 ^ 2 * t1712 + t1700 ^ 2 * t1714) * t1909 + (t1703 ^ 2 * t1710 + t1704 ^ 2 * t1712 + t1705 ^ 2 * t1714) * t1907 + ((t1784 * t1817 + t1785 * t1816 + t1786 * t1815) * t1858 + (t1787 * t1817 + t1788 * t1816 + t1789 * t1815) * t1857 + (t1675 * t1902 + t1676 * t1904 + t1678 * t1903) * MDP(10) + (t1677 * t1904 + t1679 * t1903 + t1681 * t1902) * MDP(11) + (t1698 * t1855 + t1699 * t1853 + t1700 * t1851) * t1933 + (-t1698 * t1856 - t1699 * t1854 - t1700 * t1852) * t1932) * t1796 + t1931 * (t1866 + t1867 + t1868); t1665; t1664; t1665; 0.3e1 * MDP(4) + MDP(12);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
