% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G3A0
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
%   see P3RPRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G3A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:12
% EndTime: 2020-03-09 21:27:15
% DurationCPUTime: 3.20s
% Computational Cost: add. (34320->249), mult. (26538->386), div. (1473->10), fcn. (13608->59), ass. (0->199)
t1873 = pkin(7) + qJ(3,3);
t1881 = sin(qJ(3,3));
t1799 = pkin(1) * sin(t1873) + pkin(2) * t1881;
t1791 = 0.1e1 / t1799 ^ 2;
t1878 = legFrame(3,2);
t1946 = qJ(1,3) + pkin(7);
t1835 = t1878 + t1946;
t1829 = qJ(3,3) + t1835;
t1836 = -t1878 + t1946;
t1830 = qJ(3,3) + t1836;
t1772 = cos(t1830) + cos(t1829);
t1841 = sin(qJ(1,3) + t1873);
t1894 = xDP(2);
t1895 = xDP(1);
t1945 = sin(t1829) - sin(t1830);
t1893 = xDP(3);
t1993 = -2 * t1893;
t1754 = t1772 * t1895 + t1841 * t1993 - t1945 * t1894;
t1969 = t1754 / 0.2e1;
t1996 = t1791 * t1969;
t1874 = pkin(7) + qJ(3,2);
t1883 = sin(qJ(3,2));
t1800 = pkin(1) * sin(t1874) + pkin(2) * t1883;
t1794 = 0.1e1 / t1800 ^ 2;
t1879 = legFrame(2,2);
t1947 = qJ(1,2) + pkin(7);
t1837 = t1879 + t1947;
t1831 = qJ(3,2) + t1837;
t1838 = -t1879 + t1947;
t1832 = qJ(3,2) + t1838;
t1773 = cos(t1832) + cos(t1831);
t1842 = sin(qJ(1,2) + t1874);
t1944 = sin(t1831) - sin(t1832);
t1756 = t1773 * t1895 + t1842 * t1993 - t1944 * t1894;
t1965 = t1756 / 0.2e1;
t1998 = t1794 * t1965;
t1875 = pkin(7) + qJ(3,1);
t1885 = sin(qJ(3,1));
t1801 = pkin(1) * sin(t1875) + pkin(2) * t1885;
t1797 = 0.1e1 / t1801 ^ 2;
t1880 = legFrame(1,2);
t1948 = qJ(1,1) + pkin(7);
t1839 = t1880 + t1948;
t1833 = qJ(3,1) + t1839;
t1840 = -t1880 + t1948;
t1834 = qJ(3,1) + t1840;
t1774 = cos(t1834) + cos(t1833);
t1843 = sin(qJ(1,1) + t1875);
t1943 = sin(t1833) - sin(t1834);
t1755 = t1774 * t1895 + t1843 * t1993 - t1943 * t1894;
t1967 = t1755 / 0.2e1;
t1997 = t1797 * t1967;
t1898 = pkin(1) ^ 2;
t1994 = MDP(4) * t1898 + MDP(1);
t2008 = t1994 / 0.2e1;
t2007 = MDP(5) / 0.2e1;
t2006 = MDP(6) / 0.2e1;
t2005 = MDP(7) / 0.4e1;
t1897 = 0.1e1 / pkin(3);
t2004 = t1897 / 0.2e1;
t1870 = qJ(1,1) + t1880;
t1871 = qJ(1,1) - t1880;
t1759 = t1943 * pkin(3) + (sin(t1839) - sin(t1840)) * pkin(2) + (sin(t1870) - sin(t1871)) * pkin(1);
t1762 = -t1774 * pkin(3) + (-cos(t1839) - cos(t1840)) * pkin(2) + (-cos(t1870) - cos(t1871)) * pkin(1);
t1768 = sin(qJ(1,1)) * pkin(1) + pkin(2) * sin(t1948) + pkin(3) * t1843;
t1796 = 0.1e1 / t1801;
t1940 = 2 * t1893;
t1975 = (t1759 * t1894 + t1762 * t1895 + t1768 * t1940) * t1796 * t1897;
t1743 = t1975 / 0.2e1;
t1749 = t1796 * t1967;
t1740 = t1749 + t1743;
t1877 = cos(pkin(7));
t1856 = pkin(1) * t1877 + pkin(2);
t1889 = cos(qJ(3,1));
t1876 = sin(pkin(7));
t1950 = t1876 * t1885;
t1783 = -pkin(1) * t1950 + t1889 * t1856;
t1992 = pkin(1) * t1876;
t1786 = t1856 * t1885 + t1889 * t1992;
t1910 = (pkin(3) + t1783) * t1740 / t1786 * t1975;
t2003 = -t1910 / 0.2e1;
t1868 = qJ(1,2) + t1879;
t1869 = qJ(1,2) - t1879;
t1758 = t1944 * pkin(3) + (sin(t1837) - sin(t1838)) * pkin(2) + (sin(t1868) - sin(t1869)) * pkin(1);
t1761 = -t1773 * pkin(3) + (-cos(t1837) - cos(t1838)) * pkin(2) + (-cos(t1868) - cos(t1869)) * pkin(1);
t1767 = sin(qJ(1,2)) * pkin(1) + pkin(2) * sin(t1947) + pkin(3) * t1842;
t1793 = 0.1e1 / t1800;
t1974 = (t1758 * t1894 + t1761 * t1895 + t1767 * t1940) * t1793 * t1897;
t1744 = t1974 / 0.2e1;
t1750 = t1793 * t1965;
t1741 = t1750 + t1744;
t1888 = cos(qJ(3,2));
t1951 = t1876 * t1883;
t1782 = -pkin(1) * t1951 + t1888 * t1856;
t1785 = t1856 * t1883 + t1888 * t1992;
t1911 = (pkin(3) + t1782) * t1741 / t1785 * t1974;
t2002 = -t1911 / 0.2e1;
t1866 = qJ(1,3) + t1878;
t1867 = qJ(1,3) - t1878;
t1757 = t1945 * pkin(3) + (sin(t1835) - sin(t1836)) * pkin(2) + (sin(t1866) - sin(t1867)) * pkin(1);
t1760 = -t1772 * pkin(3) + (-cos(t1835) - cos(t1836)) * pkin(2) + (-cos(t1866) - cos(t1867)) * pkin(1);
t1766 = sin(qJ(1,3)) * pkin(1) + pkin(2) * sin(t1946) + pkin(3) * t1841;
t1790 = 0.1e1 / t1799;
t1976 = (t1757 * t1894 + t1760 * t1895 + t1766 * t1940) * t1790 * t1897;
t1742 = t1976 / 0.2e1;
t1748 = t1790 * t1969;
t1739 = t1748 + t1742;
t1887 = cos(qJ(3,3));
t1952 = t1876 * t1881;
t1781 = -pkin(1) * t1952 + t1887 * t1856;
t1784 = t1856 * t1881 + t1887 * t1992;
t1912 = (pkin(3) + t1781) * t1739 / t1784 * t1976;
t2001 = -t1912 / 0.2e1;
t1941 = 0.2e1 * MDP(7);
t2000 = t1941 / 0.2e1;
t1942 = 0.2e1 * MDP(6);
t1999 = t1942 / 0.4e1;
t1995 = pkin(2) * t1877;
t1991 = pkin(3) * t1739;
t1990 = pkin(3) * t1741;
t1989 = t1740 * pkin(3);
t1988 = t1887 * pkin(2);
t1987 = t1888 * pkin(2);
t1986 = t1889 * pkin(2);
t1985 = 0.2e1 * pkin(2) * pkin(3);
t1983 = 0.2e1 * pkin(1);
t1729 = t1793 * t1744 * t1990;
t1864 = cos(t1874);
t1730 = -t1990 + (-pkin(1) * t1864 - t1987) * t1750;
t1735 = t1750 + t1744 / 0.2e1;
t1872 = pkin(2) ^ 2 + t1898;
t1896 = pkin(3) ^ 2;
t1977 = (t1735 * t1888 * t1985 + t1872 * t1750 + t1741 * t1896 + (pkin(3) * t1735 * t1864 + t1750 * t1995) * t1983) * t1897;
t1715 = t1729 + t2002 + (-t1730 - t1977) * t1998;
t1982 = t1715 * t1793;
t1727 = t1796 * t1743 * t1989;
t1865 = cos(t1875);
t1731 = -t1989 + (-pkin(1) * t1865 - t1986) * t1749;
t1734 = t1749 + t1743 / 0.2e1;
t1978 = (t1734 * t1889 * t1985 + t1872 * t1749 + t1740 * t1896 + (pkin(3) * t1734 * t1865 + t1749 * t1995) * t1983) * t1897;
t1717 = t1727 + t2003 + (-t1731 - t1978) * t1997;
t1981 = t1717 * t1796;
t1728 = t1790 * t1742 * t1991;
t1863 = cos(t1873);
t1732 = -t1991 + (-pkin(1) * t1863 - t1988) * t1748;
t1733 = t1748 + t1742 / 0.2e1;
t1979 = (t1733 * t1887 * t1985 + t1872 * t1748 + t1739 * t1896 + (pkin(3) * t1733 * t1863 + t1748 * t1995) * t1983) * t1897;
t1719 = t1728 + t2001 + (-t1732 - t1979) * t1996;
t1980 = t1719 * t1790;
t1973 = t1754 ^ 2 * t1790 * t1791;
t1972 = t1755 ^ 2 * t1796 * t1797;
t1971 = t1756 ^ 2 * t1793 * t1794;
t1964 = t1945 * t1790;
t1963 = t1944 * t1793;
t1962 = t1943 * t1796;
t1961 = t1772 * t1790;
t1960 = t1773 * t1793;
t1959 = t1774 * t1796;
t1958 = t1784 * t1790;
t1957 = t1785 * t1793;
t1956 = t1786 * t1796;
t1955 = t1790 * t1841;
t1954 = t1793 * t1842;
t1953 = t1796 * t1843;
t1949 = t1897 / 0.8e1;
t1939 = (t1729 - t1911 / 0.4e1 + (-t1730 - t1977 / 0.2e1) * t1998) * t1957;
t1938 = (t1727 - t1910 / 0.4e1 + (-t1731 - t1978 / 0.2e1) * t1997) * t1956;
t1937 = (t1728 - t1912 / 0.4e1 + (-t1732 - t1979 / 0.2e1) * t1996) * t1958;
t1763 = -t1987 + (-t1877 * t1888 + t1951) * pkin(1);
t1936 = (0.2e1 * t1729 + t2002 + (-0.2e1 * t1730 - t1977) * t1998) * t1763 * t1793;
t1765 = t1986 + (t1877 * t1889 - t1950) * pkin(1);
t1935 = (0.2e1 * t1727 + t2003 + (-0.2e1 * t1731 - t1978) * t1997) * t1765 * t1796;
t1764 = t1988 + (t1877 * t1887 - t1952) * pkin(1);
t1934 = (0.2e1 * t1728 + t2001 + (-0.2e1 * t1732 - t1979) * t1996) * t1764 * t1790;
t1724 = -t1730 * t1998 + t1729;
t1933 = t1724 * t1782 * t1793;
t1932 = t1724 * t1957;
t1725 = -t1731 * t1997 + t1727;
t1931 = t1725 * t1783 * t1796;
t1930 = t1725 * t1956;
t1726 = -t1732 * t1996 + t1728;
t1929 = t1726 * t1781 * t1790;
t1928 = t1726 * t1958;
t1927 = t1790 * t1976;
t1926 = t1796 * t1975;
t1925 = t1793 * t1974;
t1924 = t1781 * t1973;
t1923 = t1784 * t1973;
t1922 = t1783 * t1972;
t1921 = t1786 * t1972;
t1920 = t1782 * t1971;
t1919 = t1785 * t1971;
t1918 = t1733 * t1784 * t1927;
t1917 = t1734 * t1786 * t1926;
t1916 = t1735 * t1785 * t1925;
t1915 = (t1742 + 0.2e1 * t1748) * t1764 * t1927;
t1914 = (t1743 + 0.2e1 * t1749) * t1765 * t1926;
t1913 = (t1744 + 0.2e1 * t1750) * t1763 * t1925;
t1 = [((t1760 * t1923 + t1761 * t1919 + t1762 * t1921) * MDP(6) + (t1760 * t1924 + t1761 * t1920 + t1762 * t1922) * MDP(7)) * t1949 + (-t1772 * t1918 - t1773 * t1916 - t1774 * t1917) * t1999 + (-t1772 * t1915 + t1773 * t1913 - t1774 * t1914) * t2005 + (t1715 * t1960 + t1717 * t1959 + t1719 * t1961) * t2007 + (t1772 * t1934 - t1773 * t1936 + t1774 * t1935) * t2006 + (-t1772 * t1937 - t1773 * t1939 - t1774 * t1938) * t2000 + ((t1760 * t1980 + t1761 * t1982 + t1762 * t1981) * MDP(5) + (t1760 * t1929 + t1761 * t1933 + t1762 * t1931) * MDP(6) + (-t1760 * t1928 - t1761 * t1932 - t1762 * t1930) * MDP(7)) * t2004 + (t1724 * t1960 + t1725 * t1959 + t1726 * t1961) * t2008; ((t1757 * t1923 + t1758 * t1919 + t1759 * t1921) * MDP(6) + (t1757 * t1924 + t1758 * t1920 + t1759 * t1922) * MDP(7)) * t1949 + (t1916 * t1944 + t1917 * t1943 + t1918 * t1945) * t1999 + (-t1913 * t1944 + t1914 * t1943 + t1915 * t1945) * t2005 + (-t1715 * t1963 - t1717 * t1962 - t1719 * t1964) * t2007 + (-t1934 * t1945 - t1935 * t1943 + t1936 * t1944) * t2006 + (t1937 * t1945 + t1938 * t1943 + t1939 * t1944) * t2000 + ((t1757 * t1980 + t1758 * t1982 + t1759 * t1981) * MDP(5) + (t1757 * t1929 + t1758 * t1933 + t1759 * t1931) * MDP(6) + (-t1757 * t1928 - t1758 * t1932 - t1759 * t1930) * MDP(7)) * t2004 + (-t1724 * t1963 - t1725 * t1962 - t1726 * t1964) * t2008; (-t1715 * t1954 - t1717 * t1953 - t1719 * t1955) * MDP(5) + (-t1841 * t1934 + t1842 * t1936 - t1843 * t1935) * MDP(6) + (t1841 * t1937 + t1842 * t1939 + t1843 * t1938) * t1941 + (t1841 * t1918 + t1842 * t1916 + t1843 * t1917) * t1942 / 0.2e1 + (t1841 * t1915 - t1842 * t1913 + t1843 * t1914) * MDP(7) / 0.2e1 + ((t1766 * t1980 + t1767 * t1982 + t1768 * t1981) * MDP(5) + (t1766 * t1929 + t1767 * t1933 + t1768 * t1931) * MDP(6) + (-t1766 * t1928 - t1767 * t1932 - t1768 * t1930) * MDP(7) + (t1766 * t1923 + t1767 * t1919 + t1768 * t1921) * MDP(6) / 0.4e1 + (t1766 * t1924 + t1767 * t1920 + t1768 * t1922) * t2005) * t1897 + t1994 * (-t1724 * t1954 - t1725 * t1953 - t1726 * t1955);];
taucX  = t1;
