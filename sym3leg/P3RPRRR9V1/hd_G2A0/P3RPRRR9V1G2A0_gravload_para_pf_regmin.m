% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR9V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:52:28
% EndTime: 2020-08-06 18:52:29
% DurationCPUTime: 1.32s
% Computational Cost: add. (852->175), mult. (1074->310), div. (153->10), fcn. (1095->26), ass. (0->161)
t1991 = pkin(5) + pkin(6);
t1877 = qJ(2,3) + t1991;
t1865 = 0.1e1 / t1877;
t1899 = cos(qJ(1,3));
t1950 = t1899 * t1865;
t1893 = sin(qJ(1,3));
t1889 = legFrame(3,2);
t1868 = sin(t1889);
t1871 = cos(t1889);
t1921 = t1871 * g(1) - t1868 * g(2);
t1996 = g(3) * t1893 - t1921 * t1899;
t1998 = t1996 * t1950;
t1878 = qJ(2,2) + t1991;
t1866 = 0.1e1 / t1878;
t1901 = cos(qJ(1,2));
t1948 = t1901 * t1866;
t1895 = sin(qJ(1,2));
t1890 = legFrame(2,2);
t1869 = sin(t1890);
t1872 = cos(t1890);
t1920 = t1872 * g(1) - t1869 * g(2);
t1995 = g(3) * t1895 - t1920 * t1901;
t2000 = t1995 * t1948;
t1879 = qJ(2,1) + t1991;
t1867 = 0.1e1 / t1879;
t1903 = cos(qJ(1,1));
t1946 = t1903 * t1867;
t1897 = sin(qJ(1,1));
t1891 = legFrame(1,2);
t1870 = sin(t1891);
t1873 = cos(t1891);
t1919 = t1873 * g(1) - t1870 * g(2);
t1994 = g(3) * t1897 - t1919 * t1903;
t2002 = t1994 * t1946;
t1912 = t2002 + t2000 + t1998;
t1881 = pkin(7) + qJ(3,3);
t1858 = sin(t1881);
t1861 = cos(t1881);
t1960 = t1868 * t1893;
t1815 = t1871 * t1858 - t1861 * t1960;
t1882 = pkin(7) + qJ(3,2);
t1859 = sin(t1882);
t1862 = cos(t1882);
t1959 = t1869 * t1895;
t1817 = t1872 * t1859 - t1862 * t1959;
t1883 = pkin(7) + qJ(3,1);
t1860 = sin(t1883);
t1863 = cos(t1883);
t1958 = t1870 * t1897;
t1819 = t1873 * t1860 - t1863 * t1958;
t1855 = 0.1e1 / t1861;
t1969 = t1855 * t1865;
t1997 = t1996 * t1969;
t1856 = 0.1e1 / t1862;
t1966 = t1856 * t1866;
t1999 = t1995 * t1966;
t1857 = 0.1e1 / t1863;
t1963 = t1857 * t1867;
t2001 = t1994 * t1963;
t1910 = t1815 * t1997 + t1817 * t1999 + t1819 * t2001;
t1957 = t1871 * t1893;
t1816 = t1868 * t1858 + t1861 * t1957;
t1956 = t1872 * t1895;
t1818 = t1869 * t1859 + t1862 * t1956;
t1955 = t1873 * t1897;
t1820 = t1870 * t1860 + t1863 * t1955;
t1908 = t1816 * t1997 + t1818 * t1999 + t1820 * t2001;
t1993 = 2 * pkin(3);
t1888 = cos(pkin(7));
t1992 = 0.2e1 * t1888 ^ 2;
t1898 = cos(qJ(3,3));
t1884 = t1898 ^ 2;
t1987 = t1884 * pkin(3);
t1900 = cos(qJ(3,2));
t1885 = t1900 ^ 2;
t1986 = t1885 * pkin(3);
t1902 = cos(qJ(3,1));
t1886 = t1902 ^ 2;
t1985 = t1886 * pkin(3);
t1984 = t1898 * pkin(2);
t1983 = t1900 * pkin(2);
t1982 = t1902 * pkin(2);
t1874 = t1897 * pkin(1);
t1806 = -g(3) * (t1903 * qJ(2,1) - t1874) - t1919 * (t1903 * pkin(1) + t1897 * qJ(2,1));
t1981 = t1806 * t1857;
t1875 = pkin(1) * t1893;
t1807 = -g(3) * (t1899 * qJ(2,3) - t1875) - t1921 * (t1899 * pkin(1) + t1893 * qJ(2,3));
t1980 = t1807 * t1855;
t1876 = pkin(1) * t1895;
t1808 = -g(3) * (t1901 * qJ(2,2) - t1876) - t1920 * (t1901 * pkin(1) + t1895 * qJ(2,2));
t1979 = t1808 * t1856;
t1978 = t1996 * t1865;
t1977 = t1995 * t1866;
t1976 = t1994 * t1867;
t1887 = sin(pkin(7));
t1892 = sin(qJ(3,3));
t1954 = t1887 * t1892;
t1975 = t1996 / (t1888 * t1898 - t1954);
t1894 = sin(qJ(3,2));
t1953 = t1887 * t1894;
t1974 = t1995 / (t1888 * t1900 - t1953);
t1896 = sin(qJ(3,1));
t1952 = t1896 * t1887;
t1973 = t1994 / (t1902 * t1888 - t1952);
t1905 = pkin(2) / 0.2e1;
t1972 = (t1898 * pkin(3) + t1905) * t1892;
t1971 = (t1900 * pkin(3) + t1905) * t1894;
t1970 = (t1902 * pkin(3) + t1905) * t1896;
t1968 = t1855 * t1868;
t1967 = t1855 * t1871;
t1965 = t1856 * t1869;
t1964 = t1856 * t1872;
t1962 = t1857 * t1870;
t1961 = t1857 * t1873;
t1864 = pkin(1) * t1887;
t1951 = t1898 * (-t1892 * pkin(3) + t1864);
t1949 = t1900 * (-t1894 * pkin(3) + t1864);
t1947 = t1902 * (-t1896 * pkin(3) + t1864);
t1915 = g(3) * t1899 + t1921 * t1893;
t1936 = t1915 * t1969;
t1914 = g(3) * t1901 + t1920 * t1895;
t1935 = t1914 * t1966;
t1913 = g(3) * t1903 + t1919 * t1897;
t1934 = t1913 * t1963;
t1933 = t1893 * t1954;
t1932 = t1895 * t1953;
t1931 = t1897 * t1952;
t1927 = -t1903 * t1879 + t1874;
t1926 = -t1899 * t1877 + t1875;
t1925 = -t1901 * t1878 + t1876;
t1924 = t1858 * t1997;
t1923 = t1859 * t1999;
t1922 = t1860 * t2001;
t1918 = pkin(2) * t1933 + (t1933 * t1993 - t1926) * t1898;
t1917 = pkin(2) * t1932 + (t1932 * t1993 - t1925) * t1900;
t1916 = pkin(2) * t1931 + (t1931 * t1993 - t1927) * t1902;
t1911 = t1913 * t1946 + t1914 * t1948 + t1915 * t1950;
t1909 = t1815 * t1936 + t1817 * t1935 + t1819 * t1934;
t1907 = t1816 * t1936 + t1818 * t1935 + t1820 * t1934;
t1906 = 1 / pkin(3);
t1904 = -pkin(3) / 0.2e1;
t1854 = t1888 * pkin(2) + pkin(1);
t1847 = t1985 + t1982 / 0.2e1 + t1904;
t1846 = t1986 + t1983 / 0.2e1 + t1904;
t1845 = t1987 + t1984 / 0.2e1 + t1904;
t1844 = t1870 * g(1) + t1873 * g(2);
t1843 = t1869 * g(1) + t1872 * g(2);
t1842 = t1868 * g(1) + t1871 * g(2);
t1814 = pkin(1) * t1896 + (-pkin(3) + t1982 + 0.2e1 * t1985) * t1887;
t1813 = pkin(1) * t1894 + (-pkin(3) + t1983 + 0.2e1 * t1986) * t1887;
t1812 = pkin(1) * t1892 + (-pkin(3) + t1984 + 0.2e1 * t1987) * t1887;
t1811 = t1927 * t1952 + (t1886 - 0.1e1) * t1897 * pkin(3);
t1810 = t1925 * t1953 + (t1885 - 0.1e1) * t1895 * pkin(3);
t1809 = t1926 * t1954 + (t1884 - 0.1e1) * t1893 * pkin(3);
t1805 = t1842 * t1858 + t1861 * t1915;
t1804 = -t1842 * t1861 + t1858 * t1915;
t1803 = t1844 * t1860 + t1863 * t1913;
t1802 = -t1844 * t1863 + t1860 * t1913;
t1801 = t1843 * t1859 + t1862 * t1914;
t1800 = -t1843 * t1862 + t1859 * t1914;
t1 = [0, t1908, t1907, t1908 * t1888, -t1908 * t1887, -t1907, (t1820 * t1981 - ((t1847 * t1955 + t1870 * t1970) * t1992 + (t1870 * t1814 - t1916 * t1873) * t1888 - t1811 * t1873 + t1870 * t1947) * t1973) * t1867 + (t1818 * t1979 - ((t1846 * t1956 + t1869 * t1971) * t1992 + (t1869 * t1813 - t1917 * t1872) * t1888 - t1810 * t1872 + t1869 * t1949) * t1974) * t1866 + (t1816 * t1980 - ((t1845 * t1957 + t1868 * t1972) * t1992 + (t1868 * t1812 - t1918 * t1871) * t1888 - t1809 * t1871 + t1868 * t1951) * t1975) * t1865, 0, 0, 0, 0, 0, t1816 * t1978 + t1818 * t1977 + t1820 * t1976 + (t1800 * t1965 + t1802 * t1962 + t1804 * t1968) * t1906, -t1816 * t1924 - t1818 * t1923 - t1820 * t1922 + (t1801 * t1965 + t1803 * t1962 + t1805 * t1968) * t1906, -g(1); 0, t1910, t1909, t1910 * t1888, -t1910 * t1887, -t1909, (t1819 * t1981 - ((-t1847 * t1958 + t1873 * t1970) * t1992 + (t1873 * t1814 + t1916 * t1870) * t1888 + t1811 * t1870 + t1873 * t1947) * t1973) * t1867 + (t1817 * t1979 - ((-t1846 * t1959 + t1872 * t1971) * t1992 + (t1872 * t1813 + t1917 * t1869) * t1888 + t1810 * t1869 + t1872 * t1949) * t1974) * t1866 + (t1815 * t1980 - ((-t1845 * t1960 + t1871 * t1972) * t1992 + (t1871 * t1812 + t1918 * t1868) * t1888 + t1809 * t1868 + t1871 * t1951) * t1975) * t1865, 0, 0, 0, 0, 0, t1815 * t1978 + t1817 * t1977 + t1819 * t1976 + (t1800 * t1964 + t1802 * t1961 + t1804 * t1967) * t1906, -t1815 * t1924 - t1817 * t1923 - t1819 * t1922 + (t1801 * t1964 + t1803 * t1961 + t1805 * t1967) * t1906, -g(2); 0, t1912, t1911, t1912 * t1888, -t1912 * t1887, -t1911, (-t1897 * t1879 * t1994 + (t1806 - (pkin(3) * t1863 + t1854) * t1994) * t1903) * t1867 + (-t1895 * t1878 * t1995 + (t1808 - (pkin(3) * t1862 + t1854) * t1995) * t1901) * t1866 + (-t1893 * t1877 * t1996 + (t1807 - (pkin(3) * t1861 + t1854) * t1996) * t1899) * t1865, 0, 0, 0, 0, 0, t1861 * t1998 + t1862 * t2000 + t1863 * t2002, -t1858 * t1998 - t1859 * t2000 - t1860 * t2002, -g(3);];
tau_reg  = t1;
