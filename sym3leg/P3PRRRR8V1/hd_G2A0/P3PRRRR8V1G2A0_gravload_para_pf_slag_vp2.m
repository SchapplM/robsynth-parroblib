% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:55
% EndTime: 2020-08-06 17:02:56
% DurationCPUTime: 0.78s
% Computational Cost: add. (675->134), mult. (1503->253), div. (54->7), fcn. (1311->22), ass. (0->117)
t880 = sin(qJ(2,1));
t874 = legFrame(1,2);
t859 = sin(t874);
t862 = cos(t874);
t842 = t862 * g(1) - t859 * g(2);
t869 = cos(pkin(6));
t924 = t869 * t842;
t867 = sin(pkin(6));
t946 = t867 * g(3);
t836 = t924 - t946;
t839 = t859 * g(1) + t862 * g(2);
t868 = sin(pkin(3));
t870 = cos(pkin(3));
t897 = t836 * t870 + t839 * t868;
t959 = t897 * t880;
t878 = sin(qJ(2,2));
t873 = legFrame(2,2);
t858 = sin(t873);
t861 = cos(t873);
t841 = t861 * g(1) - t858 * g(2);
t925 = t869 * t841;
t835 = t925 - t946;
t838 = t858 * g(1) + t861 * g(2);
t898 = t835 * t870 + t838 * t868;
t958 = t898 * t878;
t876 = sin(qJ(2,3));
t872 = legFrame(3,2);
t857 = sin(t872);
t860 = cos(t872);
t840 = t860 * g(1) - t857 * g(2);
t926 = t869 * t840;
t834 = t926 - t946;
t837 = t857 * g(1) + t860 * g(2);
t899 = t834 * t870 + t837 * t868;
t957 = t899 * t876;
t882 = cos(qJ(2,3));
t881 = cos(qJ(3,3));
t950 = pkin(2) * t881;
t843 = -t882 * pkin(5) + t876 * t950;
t875 = sin(qJ(3,3));
t953 = pkin(2) * t875;
t831 = t843 * t868 + t870 * t953;
t956 = 0.1e1 / t831;
t884 = cos(qJ(2,2));
t883 = cos(qJ(3,2));
t949 = pkin(2) * t883;
t844 = -t884 * pkin(5) + t878 * t949;
t877 = sin(qJ(3,2));
t952 = pkin(2) * t877;
t832 = t844 * t868 + t870 * t952;
t955 = 0.1e1 / t832;
t886 = cos(qJ(2,1));
t885 = cos(qJ(3,1));
t948 = pkin(2) * t885;
t845 = -t886 * pkin(5) + t880 * t948;
t879 = sin(qJ(3,1));
t951 = pkin(2) * t879;
t833 = t845 * t868 + t870 * t951;
t954 = 0.1e1 / t833;
t947 = g(3) * t869;
t945 = t956 / t881;
t944 = t955 / t883;
t943 = t954 / t885;
t942 = t956 * t837;
t941 = t955 * t838;
t940 = t954 * t839;
t935 = t837 * t870;
t933 = t838 * t870;
t931 = t839 * t870;
t930 = mrSges(3,2) * t946 * t868;
t929 = t867 * t882;
t928 = t867 * t884;
t927 = t867 * t886;
t923 = t870 * t876;
t922 = t870 * t878;
t921 = t870 * t880;
t920 = t870 * t882;
t919 = t870 * t884;
t918 = t870 * t886;
t917 = t875 * t882;
t916 = t877 * t884;
t915 = t879 * t886;
t902 = t840 * t867 + t947;
t890 = t902 * t882 + t957;
t914 = (((t834 * t868 - t935) * mrSges(3,1) + t890 * mrSges(3,2)) * t881 + (t930 + (-t868 * t926 + t935) * mrSges(3,2) + t890 * mrSges(3,1)) * t875) * t945;
t901 = t841 * t867 + t947;
t889 = t901 * t884 + t958;
t913 = (((t835 * t868 - t933) * mrSges(3,1) + t889 * mrSges(3,2)) * t883 + (t930 + (-t868 * t925 + t933) * mrSges(3,2) + t889 * mrSges(3,1)) * t877) * t944;
t900 = t842 * t867 + t947;
t888 = t900 * t886 + t959;
t912 = (((t836 * t868 - t931) * mrSges(3,1) + t888 * mrSges(3,2)) * t885 + (t930 + (-t868 * t924 + t931) * mrSges(3,2) + t888 * mrSges(3,1)) * t879) * t943;
t871 = mrSges(2,2) - mrSges(3,3);
t852 = t871 * t947;
t911 = (t852 * t882 + (t840 * t929 + t957) * t871 + (t902 * t876 - t899 * t882) * (t881 * mrSges(3,1) - t875 * mrSges(3,2) + mrSges(2,1))) * t945;
t910 = (t852 * t884 + (t841 * t928 + t958) * t871 + (t901 * t878 - t898 * t884) * (t883 * mrSges(3,1) - t877 * mrSges(3,2) + mrSges(2,1))) * t944;
t909 = (t852 * t886 + (t842 * t927 + t959) * t871 + (t900 * t880 - t897 * t886) * (t885 * mrSges(3,1) - t879 * mrSges(3,2) + mrSges(2,1))) * t943;
t908 = ((t867 * t920 + t869 * t876) * t950 + (t867 * t923 - t869 * t882) * pkin(5)) * t914;
t907 = ((t867 * t919 + t869 * t878) * t949 + (t867 * t922 - t869 * t884) * pkin(5)) * t913;
t906 = ((t867 * t918 + t869 * t880) * t948 + (t867 * t921 - t869 * t886) * pkin(5)) * t912;
t893 = t868 * t881 + t875 * t923;
t905 = (t893 * t867 - t869 * t917) * t911;
t892 = t868 * t883 + t877 * t922;
t904 = (t892 * t867 - t869 * t916) * t910;
t891 = t868 * t885 + t879 * t921;
t903 = (t891 * t867 - t869 * t915) * t909;
t896 = -t843 * t870 + t868 * t953;
t895 = -t844 * t870 + t868 * t952;
t894 = -t845 * t870 + t868 * t951;
t887 = 0.1e1 / pkin(2);
t863 = m(1) + m(2) + m(3);
t848 = pkin(5) * t880 + t886 * t948;
t847 = pkin(5) * t878 + t884 * t949;
t846 = pkin(5) * t876 + t882 * t950;
t818 = -t867 * t848 + t894 * t869;
t817 = -t867 * t847 + t895 * t869;
t816 = -t867 * t846 + t896 * t869;
t1 = [-t860 * t905 - t861 * t904 - t862 * t903 - g(1) * m(4) + (-t860 * t908 - t861 * t907 - t862 * t906) * t887 + (-(-t818 * t862 + t833 * t859) * t940 - (-t817 * t861 + t832 * t858) * t941 - (-t816 * t860 + t831 * t857) * t942) * t863; t857 * t905 + t858 * t904 + t859 * t903 - g(2) * m(4) + (t857 * t908 + t858 * t907 + t859 * t906) * t887 + (-(t818 * t859 + t833 * t862) * t940 - (t817 * t858 + t832 * t861) * t941 - (t816 * t857 + t831 * t860) * t942) * t863; (-t867 * t915 - t891 * t869) * t909 + (-t867 * t916 - t892 * t869) * t910 + (-t867 * t917 - t893 * t869) * t911 - g(3) * m(4) + ((-(-t867 * t880 + t869 * t918) * t948 - pkin(5) * (t869 * t921 + t927)) * t912 + (-(-t867 * t878 + t869 * t919) * t949 - pkin(5) * (t869 * t922 + t928)) * t913 + (-(-t867 * t876 + t869 * t920) * t950 - pkin(5) * (t869 * t923 + t929)) * t914) * t887 + (-(t869 * t848 + t894 * t867) * t940 - (t869 * t847 + t895 * t867) * t941 - (t869 * t846 + t896 * t867) * t942) * t863;];
taugX  = t1;
