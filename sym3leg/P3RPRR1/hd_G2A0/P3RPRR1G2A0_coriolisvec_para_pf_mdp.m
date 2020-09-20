% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:14
% EndTime: 2020-03-09 21:25:17
% DurationCPUTime: 3.37s
% Computational Cost: add. (34320->249), mult. (26538->385), div. (1473->10), fcn. (13608->59), ass. (0->197)
t1842 = pkin(7) + qJ(3,3);
t1850 = sin(qJ(3,3));
t1768 = pkin(1) * sin(t1842) + t1850 * pkin(2);
t1760 = 0.1e1 / t1768 ^ 2;
t1847 = legFrame(3,2);
t1912 = qJ(1,3) + pkin(7);
t1804 = t1847 + t1912;
t1798 = qJ(3,3) + t1804;
t1805 = -t1847 + t1912;
t1799 = qJ(3,3) + t1805;
t1738 = sin(t1798) + sin(t1799);
t1741 = -cos(t1799) + cos(t1798);
t1816 = cos(qJ(1,3) + t1842);
t1863 = xDP(2);
t1864 = xDP(1);
t1974 = 2 * xDP(3);
t1723 = t1738 * t1864 + t1741 * t1863 + t1816 * t1974;
t1935 = t1723 / 0.2e1;
t1962 = t1760 * t1935;
t1843 = pkin(7) + qJ(3,2);
t1851 = sin(qJ(3,2));
t1769 = pkin(1) * sin(t1843) + t1851 * pkin(2);
t1763 = 0.1e1 / t1769 ^ 2;
t1848 = legFrame(2,2);
t1913 = qJ(1,2) + pkin(7);
t1806 = t1848 + t1913;
t1800 = qJ(3,2) + t1806;
t1807 = -t1848 + t1913;
t1801 = qJ(3,2) + t1807;
t1739 = sin(t1800) + sin(t1801);
t1742 = -cos(t1801) + cos(t1800);
t1817 = cos(qJ(1,2) + t1843);
t1725 = t1739 * t1864 + t1742 * t1863 + t1817 * t1974;
t1931 = t1725 / 0.2e1;
t1964 = t1763 * t1931;
t1844 = pkin(7) + qJ(3,1);
t1852 = sin(qJ(3,1));
t1770 = pkin(1) * sin(t1844) + t1852 * pkin(2);
t1766 = 0.1e1 / t1770 ^ 2;
t1849 = legFrame(1,2);
t1914 = qJ(1,1) + pkin(7);
t1808 = t1849 + t1914;
t1802 = qJ(3,1) + t1808;
t1809 = -t1849 + t1914;
t1803 = qJ(3,1) + t1809;
t1740 = sin(t1802) + sin(t1803);
t1743 = -cos(t1803) + cos(t1802);
t1818 = cos(qJ(1,1) + t1844);
t1724 = t1740 * t1864 + t1743 * t1863 + t1818 * t1974;
t1933 = t1724 / 0.2e1;
t1963 = t1766 * t1933;
t1867 = pkin(1) ^ 2;
t1960 = MDP(4) * t1867 + MDP(1);
t1975 = t1960 / 0.2e1;
t1973 = MDP(5) / 0.2e1;
t1972 = MDP(6) / 0.2e1;
t1971 = MDP(7) / 0.4e1;
t1866 = 0.1e1 / pkin(3);
t1970 = t1866 / 0.2e1;
t1837 = qJ(1,2) + t1848;
t1838 = qJ(1,2) - t1848;
t1727 = -t1739 * pkin(3) + (-sin(t1806) - sin(t1807)) * pkin(2) + (-sin(t1837) - sin(t1838)) * pkin(1);
t1730 = -t1742 * pkin(3) + (-cos(t1806) + cos(t1807)) * pkin(2) + (-cos(t1837) + cos(t1838)) * pkin(1);
t1736 = -cos(qJ(1,2)) * pkin(1) - pkin(2) * cos(t1913) - pkin(3) * t1817;
t1762 = 0.1e1 / t1769;
t1942 = (t1727 * t1864 + t1730 * t1863 + t1736 * t1974) * t1762 * t1866;
t1711 = t1942 / 0.2e1;
t1719 = t1762 * t1931;
t1710 = t1719 + t1711;
t1846 = cos(pkin(7));
t1825 = t1846 * pkin(1) + pkin(2);
t1855 = cos(qJ(3,2));
t1845 = sin(pkin(7));
t1917 = t1845 * t1851;
t1751 = -pkin(1) * t1917 + t1855 * t1825;
t1958 = pkin(1) * t1845;
t1754 = t1851 * t1825 + t1855 * t1958;
t1879 = t1710 * (pkin(3) + t1751) / t1754 * t1942;
t1969 = -t1879 / 0.2e1;
t1839 = qJ(1,1) + t1849;
t1840 = qJ(1,1) - t1849;
t1728 = -t1740 * pkin(3) + (-sin(t1808) - sin(t1809)) * pkin(2) + (-sin(t1839) - sin(t1840)) * pkin(1);
t1731 = -t1743 * pkin(3) + (-cos(t1808) + cos(t1809)) * pkin(2) + (-cos(t1839) + cos(t1840)) * pkin(1);
t1737 = -cos(qJ(1,1)) * pkin(1) - pkin(2) * cos(t1914) - pkin(3) * t1818;
t1765 = 0.1e1 / t1770;
t1940 = (t1728 * t1864 + t1731 * t1863 + t1737 * t1974) * t1765 * t1866;
t1713 = t1940 / 0.2e1;
t1718 = t1765 * t1933;
t1709 = t1718 + t1713;
t1857 = cos(qJ(3,1));
t1916 = t1845 * t1852;
t1752 = -pkin(1) * t1916 + t1857 * t1825;
t1755 = t1852 * t1825 + t1857 * t1958;
t1880 = t1709 * (pkin(3) + t1752) / t1755 * t1940;
t1968 = -t1880 / 0.2e1;
t1835 = qJ(1,3) + t1847;
t1836 = qJ(1,3) - t1847;
t1726 = -t1738 * pkin(3) + (-sin(t1804) - sin(t1805)) * pkin(2) + (-sin(t1835) - sin(t1836)) * pkin(1);
t1729 = -t1741 * pkin(3) + (-cos(t1804) + cos(t1805)) * pkin(2) + (-cos(t1835) + cos(t1836)) * pkin(1);
t1735 = -cos(qJ(1,3)) * pkin(1) - pkin(2) * cos(t1912) - pkin(3) * t1816;
t1759 = 0.1e1 / t1768;
t1941 = (t1726 * t1864 + t1729 * t1863 + t1735 * t1974) * t1759 * t1866;
t1712 = t1941 / 0.2e1;
t1717 = t1759 * t1935;
t1708 = t1717 + t1712;
t1853 = cos(qJ(3,3));
t1918 = t1845 * t1850;
t1750 = -pkin(1) * t1918 + t1853 * t1825;
t1753 = t1850 * t1825 + t1853 * t1958;
t1881 = t1708 * (pkin(3) + t1750) / t1753 * t1941;
t1967 = -t1881 / 0.2e1;
t1910 = 0.2e1 * MDP(7);
t1966 = t1910 / 0.2e1;
t1911 = 0.2e1 * MDP(6);
t1965 = t1911 / 0.4e1;
t1961 = pkin(2) * t1846;
t1957 = pkin(2) * t1853;
t1956 = pkin(3) * t1708;
t1955 = pkin(3) * t1710;
t1954 = t1709 * pkin(3);
t1953 = t1855 * pkin(2);
t1952 = t1857 * pkin(2);
t1951 = 0.2e1 * pkin(2) * pkin(3);
t1949 = 0.2e1 * pkin(1);
t1697 = t1759 * t1712 * t1956;
t1829 = cos(t1842);
t1699 = -t1956 + (-pkin(1) * t1829 - t1957) * t1717;
t1702 = t1717 + t1712 / 0.2e1;
t1841 = pkin(2) ^ 2 + t1867;
t1865 = pkin(3) ^ 2;
t1945 = (t1702 * t1853 * t1951 + t1841 * t1717 + t1708 * t1865 + (pkin(3) * t1702 * t1829 + t1717 * t1961) * t1949) * t1866;
t1684 = t1697 + t1967 + (-t1699 - t1945) * t1962;
t1948 = t1684 * t1759;
t1698 = t1762 * t1711 * t1955;
t1830 = cos(t1843);
t1700 = -t1955 + (-pkin(1) * t1830 - t1953) * t1719;
t1704 = t1719 + t1711 / 0.2e1;
t1943 = (t1704 * t1855 * t1951 + t1841 * t1719 + t1710 * t1865 + (pkin(3) * t1704 * t1830 + t1719 * t1961) * t1949) * t1866;
t1686 = t1698 + t1969 + (-t1700 - t1943) * t1964;
t1947 = t1686 * t1762;
t1696 = t1765 * t1713 * t1954;
t1831 = cos(t1844);
t1701 = -t1954 + (-pkin(1) * t1831 - t1952) * t1718;
t1703 = t1718 + t1713 / 0.2e1;
t1944 = (t1703 * t1857 * t1951 + t1841 * t1718 + t1709 * t1865 + (pkin(3) * t1703 * t1831 + t1718 * t1961) * t1949) * t1866;
t1688 = t1696 + t1968 + (-t1701 - t1944) * t1963;
t1946 = t1688 * t1765;
t1939 = t1723 ^ 2 * t1759 * t1760;
t1938 = t1724 ^ 2 * t1765 * t1766;
t1937 = t1725 ^ 2 * t1762 * t1763;
t1930 = t1738 * t1759;
t1929 = t1739 * t1762;
t1928 = t1740 * t1765;
t1927 = t1741 * t1759;
t1926 = t1742 * t1762;
t1925 = t1743 * t1765;
t1924 = t1753 * t1759;
t1923 = t1754 * t1762;
t1922 = t1755 * t1765;
t1921 = t1759 * t1816;
t1920 = t1762 * t1817;
t1919 = t1765 * t1818;
t1915 = t1866 / 0.8e1;
t1693 = -t1699 * t1962 + t1697;
t1908 = t1693 * t1750 * t1759;
t1907 = t1693 * t1924;
t1694 = -t1700 * t1964 + t1698;
t1906 = t1694 * t1751 * t1762;
t1905 = t1694 * t1923;
t1695 = -t1701 * t1963 + t1696;
t1904 = t1695 * t1752 * t1765;
t1903 = t1695 * t1922;
t1902 = t1762 * t1942;
t1901 = t1759 * t1941;
t1900 = t1765 * t1940;
t1899 = t1750 * t1939;
t1898 = t1753 * t1939;
t1897 = t1752 * t1938;
t1896 = t1755 * t1938;
t1895 = t1751 * t1937;
t1894 = t1754 * t1937;
t1732 = -t1957 + (-t1846 * t1853 + t1918) * pkin(1);
t1893 = t1759 * (0.2e1 * t1697 + t1967 + (-0.2e1 * t1699 - t1945) * t1962) * t1732;
t1892 = (t1697 - t1881 / 0.4e1 + (-t1699 - t1945 / 0.2e1) * t1962) * t1924;
t1733 = -t1953 + (-t1846 * t1855 + t1917) * pkin(1);
t1891 = t1762 * (0.2e1 * t1698 + t1969 + (-0.2e1 * t1700 - t1943) * t1964) * t1733;
t1890 = (t1698 - t1879 / 0.4e1 + (-t1700 - t1943 / 0.2e1) * t1964) * t1923;
t1734 = t1952 + (t1846 * t1857 - t1916) * pkin(1);
t1889 = t1765 * (0.2e1 * t1696 + t1968 + (-0.2e1 * t1701 - t1944) * t1963) * t1734;
t1888 = (t1696 - t1880 / 0.4e1 + (-t1701 - t1944 / 0.2e1) * t1963) * t1922;
t1887 = t1702 * t1753 * t1901;
t1886 = t1703 * t1755 * t1900;
t1885 = t1704 * t1754 * t1902;
t1884 = (t1711 + 0.2e1 * t1719) * t1733 * t1902;
t1883 = (t1712 + 0.2e1 * t1717) * t1732 * t1901;
t1882 = (t1713 + 0.2e1 * t1718) * t1734 * t1900;
t1 = [((t1726 * t1898 + t1727 * t1894 + t1728 * t1896) * MDP(6) + (t1726 * t1899 + t1727 * t1895 + t1728 * t1897) * MDP(7)) * t1915 + (-t1738 * t1887 - t1739 * t1885 - t1740 * t1886) * t1965 + (t1738 * t1883 + t1739 * t1884 - t1740 * t1882) * t1971 + (t1684 * t1930 + t1686 * t1929 + t1688 * t1928) * t1973 + (-t1738 * t1893 - t1739 * t1891 + t1740 * t1889) * t1972 + (-t1738 * t1892 - t1739 * t1890 - t1740 * t1888) * t1966 + ((t1726 * t1948 + t1727 * t1947 + t1728 * t1946) * MDP(5) + (t1726 * t1908 + t1727 * t1906 + t1728 * t1904) * MDP(6) + (-t1726 * t1907 - t1727 * t1905 - t1728 * t1903) * MDP(7)) * t1970 + (t1693 * t1930 + t1694 * t1929 + t1695 * t1928) * t1975; ((t1729 * t1898 + t1730 * t1894 + t1731 * t1896) * MDP(6) + (t1729 * t1899 + t1730 * t1895 + t1731 * t1897) * MDP(7)) * t1915 + (-t1741 * t1887 - t1742 * t1885 - t1743 * t1886) * t1965 + (t1741 * t1883 + t1742 * t1884 - t1743 * t1882) * t1971 + (t1684 * t1927 + t1686 * t1926 + t1688 * t1925) * t1973 + (-t1741 * t1893 - t1742 * t1891 + t1743 * t1889) * t1972 + (-t1741 * t1892 - t1742 * t1890 - t1743 * t1888) * t1966 + ((t1729 * t1948 + t1730 * t1947 + t1731 * t1946) * MDP(5) + (t1729 * t1908 + t1730 * t1906 + t1731 * t1904) * MDP(6) + (-t1729 * t1907 - t1730 * t1905 - t1731 * t1903) * MDP(7)) * t1970 + (t1693 * t1927 + t1694 * t1926 + t1695 * t1925) * t1975; (t1684 * t1921 + t1686 * t1920 + t1688 * t1919) * MDP(5) + (-t1816 * t1893 - t1817 * t1891 + t1818 * t1889) * MDP(6) + (-t1816 * t1892 - t1817 * t1890 - t1818 * t1888) * t1910 + (-t1816 * t1887 - t1817 * t1885 - t1818 * t1886) * t1911 / 0.2e1 + (t1816 * t1883 + t1817 * t1884 - t1818 * t1882) * MDP(7) / 0.2e1 + ((t1735 * t1948 + t1736 * t1947 + t1737 * t1946) * MDP(5) + (t1735 * t1908 + t1736 * t1906 + t1737 * t1904) * MDP(6) + (-t1735 * t1907 - t1736 * t1905 - t1737 * t1903) * MDP(7) + (t1735 * t1898 + t1736 * t1894 + t1737 * t1896) * MDP(6) / 0.4e1 + (t1735 * t1899 + t1736 * t1895 + t1737 * t1897) * t1971) * t1866 + t1960 * (t1693 * t1921 + t1694 * t1920 + t1695 * t1919);];
taucX  = t1;
