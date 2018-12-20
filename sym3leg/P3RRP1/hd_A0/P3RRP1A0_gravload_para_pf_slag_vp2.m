% Calculate Gravitation load for parallel robot
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m [4x1]
%   mass of all robot links (including platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 18:08
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3RRP1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 18:02:21
% EndTime: 2018-12-20 18:02:24
% DurationCPUTime: 2.45s
% Computational Cost: add. (2292->347), mult. (4616->509), div. (36->3), fcn. (2066->20), ass. (0->247)
t1849 = (qJ(3,3) ^ 2);
t1859 = pkin(1) ^ 2;
t1843 = cos(qJ(1,3));
t1827 = t1843 ^ 2;
t1896 = t1827 * t1859;
t1946 = (2 * t1849) + t1859 - t1896;
t1850 = (qJ(3,2) ^ 2);
t1844 = cos(qJ(1,2));
t1828 = t1844 ^ 2;
t1895 = t1828 * t1859;
t1945 = (2 * t1850) + t1859 - t1895;
t1851 = (qJ(3,1) ^ 2);
t1845 = cos(qJ(1,1));
t1829 = t1845 ^ 2;
t1894 = t1829 * t1859;
t1944 = (2 * t1851) + t1859 - t1894;
t1943 = 2 * pkin(1);
t1837 = legFrame(3,3);
t1810 = sin(t1837);
t1788 = t1810 * qJ(3,3);
t1784 = pkin(2) * t1788;
t1942 = -0.2e1 * t1784;
t1838 = legFrame(2,3);
t1811 = sin(t1838);
t1789 = t1811 * qJ(3,2);
t1785 = pkin(2) * t1789;
t1941 = -0.2e1 * t1785;
t1839 = legFrame(1,3);
t1812 = sin(t1839);
t1790 = t1812 * qJ(3,1);
t1786 = pkin(2) * t1790;
t1940 = -0.2e1 * t1786;
t1840 = sin(qJ(1,3));
t1939 = -0.2e1 * t1840;
t1841 = sin(qJ(1,2));
t1938 = -0.2e1 * t1841;
t1842 = sin(qJ(1,1));
t1937 = -0.2e1 * t1842;
t1936 = -0.2e1 * t1843;
t1935 = -0.2e1 * t1844;
t1934 = -0.2e1 * t1845;
t1933 = pkin(2) * qJ(3,1);
t1932 = pkin(2) * qJ(3,2);
t1931 = pkin(2) * qJ(3,3);
t1930 = mrSges(3,3) - mrSges(2,2);
t1813 = cos(t1837);
t1929 = pkin(2) * t1813;
t1814 = cos(t1838);
t1928 = pkin(2) * t1814;
t1815 = cos(t1839);
t1927 = pkin(2) * t1815;
t1926 = pkin(2) * t1843;
t1925 = pkin(2) * t1844;
t1924 = pkin(2) * t1845;
t1791 = t1810 * pkin(2);
t1792 = t1811 * pkin(2);
t1793 = t1812 * pkin(2);
t1923 = qJ(3,1) * t1815;
t1922 = qJ(3,1) * t1842;
t1921 = qJ(3,2) * t1814;
t1920 = qJ(3,2) * t1841;
t1919 = qJ(3,3) * t1813;
t1918 = qJ(3,3) * t1840;
t1830 = qJ(1,3) + qJ(2,3);
t1804 = sin(t1830);
t1917 = t1804 * qJ(3,3);
t1831 = qJ(1,2) + qJ(2,2);
t1805 = sin(t1831);
t1916 = t1805 * qJ(3,2);
t1832 = qJ(1,1) + qJ(2,1);
t1806 = sin(t1832);
t1915 = t1806 * qJ(3,1);
t1807 = cos(t1830);
t1798 = t1807 ^ 2;
t1858 = pkin(2) ^ 2;
t1816 = -t1849 + t1858;
t1801 = 0.1e1 + t1816;
t1833 = 0.2e1 + t1859;
t1892 = t1840 * t1843;
t1860 = -t1801 * t1827 + 0.2e1 * t1892 * t1931;
t1887 = -0.1e1 / 0.2e1 - t1858 / 0.2e1;
t1709 = 0.1e1 / (((t1801 * t1892 + (0.2e1 * t1827 - 0.1e1) * t1931) * pkin(1) * t1804 + t1918) * t1807 * t1943 + pkin(1) * t1917 * t1936 - t1849 * t1833 + (0.2e1 * (t1849 / 0.2e1 - t1860 + t1887) * t1798 + t1860) * t1859);
t1766 = -t1810 * g(1) + t1813 * g(2);
t1769 = t1813 * g(1) + t1810 * g(2);
t1730 = -t1807 * t1766 + t1804 * t1769;
t1914 = t1709 * t1730;
t1808 = cos(t1831);
t1799 = t1808 ^ 2;
t1817 = -t1850 + t1858;
t1802 = 0.1e1 + t1817;
t1890 = t1841 * t1844;
t1861 = -t1802 * t1828 + 0.2e1 * t1890 * t1932;
t1710 = 0.1e1 / (((t1802 * t1890 + (0.2e1 * t1828 - 0.1e1) * t1932) * pkin(1) * t1805 + t1920) * t1808 * t1943 + pkin(1) * t1916 * t1935 - t1850 * t1833 + (0.2e1 * (t1850 / 0.2e1 - t1861 + t1887) * t1799 + t1861) * t1859);
t1767 = -t1811 * g(1) + t1814 * g(2);
t1770 = t1814 * g(1) + t1811 * g(2);
t1731 = -t1808 * t1767 + t1805 * t1770;
t1913 = t1710 * t1731;
t1809 = cos(t1832);
t1800 = t1809 ^ 2;
t1818 = -t1851 + t1858;
t1803 = 0.1e1 + t1818;
t1888 = t1842 * t1845;
t1862 = -t1803 * t1829 + 0.2e1 * t1888 * t1933;
t1711 = 0.1e1 / (((t1803 * t1888 + (0.2e1 * t1829 - 0.1e1) * t1933) * pkin(1) * t1806 + t1922) * t1809 * t1943 + pkin(1) * t1915 * t1934 - t1851 * t1833 + (0.2e1 * (t1851 / 0.2e1 - t1862 + t1887) * t1800 + t1862) * t1859);
t1768 = -t1812 * g(1) + t1815 * g(2);
t1771 = t1815 * g(1) + t1812 * g(2);
t1732 = -t1809 * t1768 + t1806 * t1771;
t1912 = t1711 * t1732;
t1911 = t1788 * t1807;
t1910 = t1789 * t1808;
t1909 = t1790 * t1809;
t1908 = t1801 * t1813;
t1907 = t1802 * t1814;
t1906 = t1803 * t1815;
t1905 = t1804 * t1807;
t1904 = t1805 * t1808;
t1903 = t1806 * t1809;
t1902 = t1810 * t1804;
t1901 = t1811 * t1805;
t1900 = t1812 * t1806;
t1899 = t1813 * t1807;
t1898 = t1814 * t1808;
t1897 = t1815 * t1809;
t1754 = t1788 + t1929;
t1893 = t1840 * t1754;
t1758 = t1789 + t1928;
t1891 = t1841 * t1758;
t1762 = t1790 + t1927;
t1889 = t1842 * t1762;
t1794 = m(3) * qJ(3,3) + t1930;
t1797 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t1712 = (-t1766 * t1797 - t1794 * t1769) * t1807 + (-t1766 * t1794 + t1797 * t1769) * t1804;
t1795 = m(3) * qJ(3,2) + t1930;
t1713 = (-t1767 * t1797 - t1795 * t1770) * t1808 + (-t1767 * t1795 + t1797 * t1770) * t1805;
t1796 = m(3) * qJ(3,1) + t1930;
t1714 = (-t1768 * t1797 - t1796 * t1771) * t1809 + (-t1768 * t1796 + t1797 * t1771) * t1806;
t1835 = t1858 / 0.2e1;
t1886 = t1835 - t1849 / 0.2e1;
t1885 = t1835 - t1850 / 0.2e1;
t1884 = t1835 - t1851 / 0.2e1;
t1883 = pkin(2) * t1923;
t1882 = pkin(2) * t1921;
t1881 = pkin(2) * t1919;
t1880 = t1813 * t1917;
t1879 = t1814 * t1916;
t1878 = t1815 * t1915;
t1877 = t1859 * t1892;
t1876 = t1859 * t1890;
t1875 = t1859 * t1888;
t1874 = -(2 * t1849) - t1896;
t1873 = -(2 * t1850) - t1895;
t1872 = -(2 * t1851) - t1894;
t1871 = t1813 * t1877;
t1870 = t1814 * t1876;
t1869 = t1815 * t1875;
t1755 = t1788 - t1929;
t1756 = t1791 + t1919;
t1868 = -t1755 * t1827 - t1756 * t1892;
t1867 = -t1755 * t1892 + t1756 * t1827;
t1759 = t1789 - t1928;
t1760 = t1792 + t1921;
t1866 = -t1759 * t1828 - t1760 * t1890;
t1865 = -t1759 * t1890 + t1760 * t1828;
t1763 = t1790 - t1927;
t1764 = t1793 + t1923;
t1864 = -t1763 * t1829 - t1764 * t1888;
t1863 = -t1763 * t1888 + t1764 * t1829;
t1857 = koppelP(1,1);
t1856 = koppelP(2,1);
t1855 = koppelP(3,1);
t1854 = koppelP(1,2);
t1853 = koppelP(2,2);
t1852 = koppelP(3,2);
t1848 = mrSges(4,1);
t1847 = mrSges(4,2);
t1846 = xP(3);
t1834 = 0.1e1 + t1859;
t1823 = cos(t1846);
t1822 = sin(t1846);
t1787 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1783 = -0.2e1 * t1883;
t1782 = -0.2e1 * t1882;
t1781 = -0.2e1 * t1881;
t1780 = 0.2e1 * t1786;
t1779 = 0.2e1 * t1785;
t1778 = 0.2e1 * t1784;
t1777 = -t1922 + t1924;
t1776 = -t1920 + t1925;
t1775 = -t1918 + t1926;
t1774 = t1812 * t1875;
t1773 = t1811 * t1876;
t1772 = t1810 * t1877;
t1765 = t1793 - t1923;
t1761 = t1792 - t1921;
t1757 = t1791 - t1919;
t1753 = -t1822 * t1854 + t1823 * t1857;
t1752 = -t1822 * t1853 + t1823 * t1856;
t1751 = -t1822 * t1852 + t1823 * t1855;
t1750 = -t1822 * t1857 - t1823 * t1854;
t1749 = -t1822 * t1856 - t1823 * t1853;
t1748 = -t1822 * t1855 - t1823 * t1852;
t1747 = t1812 * t1818 + t1783;
t1746 = t1811 * t1817 + t1782;
t1745 = t1810 * t1816 + t1781;
t1744 = t1812 * t1803 + t1783;
t1743 = t1811 * t1802 + t1782;
t1742 = t1810 * t1801 + t1781;
t1741 = t1744 * t1845;
t1740 = t1743 * t1844;
t1739 = t1742 * t1843;
t1738 = t1765 * t1845 - t1889;
t1737 = t1762 * t1845 + t1842 * t1765;
t1736 = t1761 * t1844 - t1891;
t1735 = t1758 * t1844 + t1841 * t1761;
t1734 = t1757 * t1843 - t1893;
t1733 = t1754 * t1843 + t1840 * t1757;
t1729 = (t1818 * t1815 + t1780) * t1845 + t1747 * t1842;
t1728 = (t1817 * t1814 + t1779) * t1844 + t1746 * t1841;
t1727 = (t1816 * t1813 + t1778) * t1843 + t1745 * t1840;
t1726 = (t1780 + t1906) * t1845 + t1842 * t1744;
t1725 = (t1779 + t1907) * t1844 + t1841 * t1743;
t1724 = (t1778 + t1908) * t1843 + t1840 * t1742;
t1723 = t1747 * t1845 + (t1884 * t1815 + t1786) * t1937;
t1722 = t1746 * t1844 + (t1885 * t1814 + t1785) * t1938;
t1721 = t1745 * t1843 + (t1886 * t1813 + t1784) * t1939;
t1708 = (t1771 * mrSges(1,2) - t1768 * t1787) * t1845 + t1842 * (t1768 * mrSges(1,2) + t1787 * t1771) + t1714;
t1707 = (t1770 * mrSges(1,2) - t1767 * t1787) * t1844 + t1841 * (t1767 * mrSges(1,2) + t1787 * t1770) + t1713;
t1706 = (t1769 * mrSges(1,2) - t1766 * t1787) * t1843 + t1840 * (t1766 * mrSges(1,2) + t1787 * t1769) + t1712;
t1705 = (-t1897 + t1900) * qJ(3,1) + (-t1741 * t1800 - t1726 * t1903 + (t1812 * t1858 + t1812 - t1883) * t1845 + (-(t1940 - t1906) * t1800 - qJ(3,1) * t1765) * t1842) * pkin(1);
t1704 = (-t1898 + t1901) * qJ(3,2) + (-t1740 * t1799 - t1725 * t1904 + (t1811 * t1858 + t1811 - t1882) * t1844 + (-(t1941 - t1907) * t1799 - qJ(3,2) * t1761) * t1841) * pkin(1);
t1703 = (-t1899 + t1902) * qJ(3,3) + (-t1739 * t1798 - t1724 * t1905 + (t1810 * t1858 + t1810 - t1881) * t1843 + (-(t1942 - t1908) * t1798 - qJ(3,3) * t1757) * t1840) * pkin(1);
t1702 = -t1878 - t1909 + (t1726 * t1800 - (t1741 + ((0.1e1 / 0.2e1 + t1884) * t1815 + t1786) * t1937) * t1903 - (t1858 * t1815 + t1786 + t1815) * t1845 + qJ(3,1) * t1889) * pkin(1);
t1701 = -t1879 - t1910 + (t1725 * t1799 - (t1740 + ((0.1e1 / 0.2e1 + t1885) * t1814 + t1785) * t1938) * t1904 - (t1858 * t1814 + t1785 + t1814) * t1844 + qJ(3,2) * t1891) * pkin(1);
t1700 = -t1880 - t1911 + (t1724 * t1798 - (t1739 + ((0.1e1 / 0.2e1 + t1886) * t1813 + t1784) * t1939) * t1905 - (t1858 * t1813 + t1784 + t1813) * t1843 + qJ(3,3) * t1893) * pkin(1);
t1699 = (-t1812 * t1944 + t1783 + t1869) * t1809 + (t1872 * t1815 + t1774 + t1780) * t1806 + (-t1738 * t1800 - t1737 * t1903 + (t1793 - 0.2e1 * t1923) * t1845 + t1842 * t1790) * pkin(1);
t1698 = (-t1811 * t1945 + t1782 + t1870) * t1808 + (t1873 * t1814 + t1773 + t1779) * t1805 + (-t1736 * t1799 - t1735 * t1904 + (t1792 - 0.2e1 * t1921) * t1844 + t1841 * t1789) * pkin(1);
t1697 = (-t1810 * t1946 + t1781 + t1871) * t1807 + (t1874 * t1813 + t1772 + t1778) * t1804 + (-t1734 * t1798 - t1733 * t1905 + (t1791 - 0.2e1 * t1919) * t1843 + t1840 * t1788) * pkin(1);
t1696 = (t1815 * t1944 + t1774 + t1940) * t1809 + (t1872 * t1812 + t1783 - t1869) * t1806 + (-t1738 * t1903 + t1737 * t1800 + t1790 * t1934 + (-t1922 - t1924) * t1815) * pkin(1);
t1695 = (t1814 * t1945 + t1773 + t1941) * t1808 + (t1873 * t1811 + t1782 - t1870) * t1805 + (-t1736 * t1904 + t1735 * t1799 + t1789 * t1935 + (-t1920 - t1925) * t1814) * pkin(1);
t1694 = (t1813 * t1946 + t1772 + t1942) * t1807 + (t1874 * t1810 + t1781 - t1871) * t1804 + (-t1734 * t1905 + t1733 * t1798 + t1788 * t1936 + (-t1918 - t1926) * t1813) * pkin(1);
t1693 = -t1834 * t1878 - t1909 + (t1723 * t1903 - t1729 * t1800 + t1777 * t1762) * pkin(1) + ((t1864 - t1927) * t1809 + t1863 * t1806) * t1859;
t1692 = -t1834 * t1879 - t1910 + (t1722 * t1904 - t1728 * t1799 + t1776 * t1758) * pkin(1) + ((t1866 - t1928) * t1808 + t1865 * t1805) * t1859;
t1691 = -t1834 * t1880 - t1911 + (t1721 * t1905 - t1727 * t1798 + t1775 * t1754) * pkin(1) + ((t1868 - t1929) * t1807 + t1867 * t1804) * t1859;
t1690 = (t1834 * t1900 - t1897) * qJ(3,1) + (t1723 * t1800 + t1729 * t1903 - t1777 * t1765) * pkin(1) + ((-t1863 + t1793) * t1809 + t1864 * t1806) * t1859;
t1689 = (t1834 * t1901 - t1898) * qJ(3,2) + (t1722 * t1799 + t1728 * t1904 - t1776 * t1761) * pkin(1) + ((-t1865 + t1792) * t1808 + t1866 * t1805) * t1859;
t1688 = (t1834 * t1902 - t1899) * qJ(3,3) + (t1721 * t1798 + t1727 * t1905 - t1775 * t1757) * pkin(1) + ((-t1867 + t1791) * t1807 + t1868 * t1804) * t1859;
t1 = [-g(1) * m(4) + (t1690 * t1714 + t1705 * t1708) * t1711 + (t1689 * t1713 + t1704 * t1707) * t1710 + (t1688 * t1712 + t1703 * t1706) * t1709 + (-t1697 * t1914 - t1698 * t1913 - t1699 * t1912) * m(3); -g(2) * m(4) + (t1693 * t1714 + t1702 * t1708) * t1711 + (t1692 * t1713 + t1701 * t1707) * t1710 + (t1691 * t1712 + t1700 * t1706) * t1709 + (-t1694 * t1914 - t1695 * t1913 - t1696 * t1912) * m(3); -(-g(1) * t1848 - g(2) * t1847) * t1822 + t1823 * (g(1) * t1847 - g(2) * t1848) + ((t1702 * t1753 + t1705 * t1750) * t1708 + (t1690 * t1750 + t1693 * t1753) * t1714 - (t1696 * t1753 + t1699 * t1750) * m(3) * t1732) * t1711 + ((t1701 * t1752 + t1704 * t1749) * t1707 + (t1689 * t1749 + t1692 * t1752) * t1713 - (t1695 * t1752 + t1698 * t1749) * m(3) * t1731) * t1710 + ((t1700 * t1751 + t1703 * t1748) * t1706 + (t1688 * t1748 + t1691 * t1751) * t1712 - (t1694 * t1751 + t1697 * t1748) * m(3) * t1730) * t1709;];
taugX  = t1;
