% Calculate Gravitation load for parallel robot
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:06
% EndTime: 2020-08-06 18:32:07
% DurationCPUTime: 0.79s
% Computational Cost: add. (770->197), mult. (875->216), div. (24->7), fcn. (462->94), ass. (0->123)
t965 = 2 * pkin(2);
t903 = legFrame(1,3);
t893 = sin(t903);
t896 = cos(t903);
t824 = -t893 * g(1) + t896 * g(2);
t827 = t896 * g(1) + t893 * g(2);
t900 = qJ(1,1) + pkin(7);
t877 = sin(t900);
t880 = cos(t900);
t964 = t824 * t877 + t827 * t880;
t902 = legFrame(2,3);
t892 = sin(t902);
t895 = cos(t902);
t823 = -t892 * g(1) + t895 * g(2);
t826 = t895 * g(1) + t892 * g(2);
t899 = qJ(1,2) + pkin(7);
t876 = sin(t899);
t879 = cos(t899);
t963 = t823 * t876 + t826 * t879;
t901 = legFrame(3,3);
t891 = sin(t901);
t894 = cos(t901);
t822 = -t891 * g(1) + t894 * g(2);
t825 = t894 * g(1) + t891 * g(2);
t898 = qJ(1,3) + pkin(7);
t875 = sin(t898);
t878 = cos(t898);
t962 = t822 * t875 + t825 * t878;
t961 = -2 * pkin(1);
t960 = -2 * pkin(2);
t911 = (-pkin(6) - pkin(5));
t958 = -2 * t911;
t957 = 2 * t911;
t912 = (m(2) + m(3));
t956 = cos(pkin(7)) * pkin(1) + pkin(2);
t955 = mrSges(3,1) * t822;
t954 = mrSges(3,1) * t823;
t953 = mrSges(3,1) * t824;
t952 = mrSges(3,2) * t822;
t951 = mrSges(3,2) * t823;
t950 = mrSges(3,2) * t824;
t949 = mrSges(3,2) * t825;
t948 = mrSges(3,2) * t826;
t947 = mrSges(3,2) * t827;
t819 = mrSges(3,1) * t825;
t864 = t912 * pkin(1) + mrSges(1,1);
t874 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3);
t924 = pkin(7) + qJ(3,3);
t884 = qJ(1,3) + t924;
t927 = -pkin(7) + qJ(3,3);
t885 = qJ(1,3) - t927;
t897 = m(3) * pkin(2) + mrSges(2,1);
t907 = cos(qJ(3,3));
t940 = 0.1e1 / (t907 * pkin(3) + t956) * ((-t949 - t955) * cos(t885) / 0.2e1 + (t819 - t952) * sin(t885) / 0.2e1 + (t949 - t955) * cos(t884) / 0.2e1 + (t819 + t952) * sin(t884) / 0.2e1 + (t822 * mrSges(1,2) + t864 * t825) * sin(qJ(1,3)) + (mrSges(1,2) * t825 - t822 * t864) * cos(qJ(1,3)) + (-t822 * t878 + t825 * t875) * t897 + t962 * t874);
t820 = mrSges(3,1) * t826;
t925 = pkin(7) + qJ(3,2);
t886 = qJ(1,2) + t925;
t928 = -pkin(7) + qJ(3,2);
t887 = qJ(1,2) - t928;
t908 = cos(qJ(3,2));
t939 = 0.1e1 / (t908 * pkin(3) + t956) * ((-t948 - t954) * cos(t887) / 0.2e1 + (t820 - t951) * sin(t887) / 0.2e1 + (t948 - t954) * cos(t886) / 0.2e1 + (t820 + t951) * sin(t886) / 0.2e1 + (t823 * mrSges(1,2) + t864 * t826) * sin(qJ(1,2)) + (mrSges(1,2) * t826 - t823 * t864) * cos(qJ(1,2)) + (-t823 * t879 + t826 * t876) * t897 + t963 * t874);
t821 = mrSges(3,1) * t827;
t926 = pkin(7) + qJ(3,1);
t888 = qJ(1,1) + t926;
t929 = -pkin(7) + qJ(3,1);
t889 = qJ(1,1) - t929;
t909 = cos(qJ(3,1));
t938 = 0.1e1 / (t909 * pkin(3) + t956) * ((-t947 - t953) * cos(t889) / 0.2e1 + (t821 - t950) * sin(t889) / 0.2e1 + (t947 - t953) * cos(t888) / 0.2e1 + (t821 + t950) * sin(t888) / 0.2e1 + (t824 * mrSges(1,2) + t864 * t827) * sin(qJ(1,1)) + (mrSges(1,2) * t827 - t824 * t864) * cos(qJ(1,1)) + (-t824 * t880 + t827 * t877) * t897 + t964 * t874);
t883 = qJ(1,1) + t903;
t882 = qJ(1,2) + t902;
t881 = qJ(1,3) + t901;
t865 = t901 + t898;
t853 = qJ(3,3) + t865;
t854 = -qJ(3,3) + t865;
t937 = sin(t853) + sin(t854);
t857 = t902 + t886;
t858 = t902 + t887;
t936 = sin(t857) + sin(t858);
t867 = t903 + t900;
t861 = qJ(3,1) + t867;
t862 = -qJ(3,1) + t867;
t935 = sin(t861) + sin(t862);
t934 = cos(t853) + cos(t854);
t933 = cos(t857) + cos(t858);
t932 = cos(t861) + cos(t862);
t931 = 2 * pkin(1);
t923 = (g(3) * t912) / 0.2e1;
t904 = sin(qJ(3,3));
t913 = 0.2e1 * qJ(3,3);
t816 = 0.1e1 / (pkin(3) * sin(t913) + t904 * t965 + (sin(t924) + sin(t927)) * pkin(1));
t916 = 0.1e1 / pkin(3);
t922 = (-g(3) * (t907 * mrSges(3,1) - t904 * mrSges(3,2)) + t962 * (mrSges(3,1) * t904 + mrSges(3,2) * t907)) * t816 * t916;
t905 = sin(qJ(3,2));
t914 = 0.2e1 * qJ(3,2);
t817 = 0.1e1 / (pkin(3) * sin(t914) + t905 * t965 + (sin(t925) + sin(t928)) * pkin(1));
t921 = (-g(3) * (t908 * mrSges(3,1) - t905 * mrSges(3,2)) + t963 * (mrSges(3,1) * t905 + mrSges(3,2) * t908)) * t817 * t916;
t906 = sin(qJ(3,1));
t915 = 0.2e1 * qJ(3,1);
t818 = 0.1e1 / (pkin(3) * sin(t915) + t906 * t965 + (sin(t926) + sin(t929)) * pkin(1));
t920 = (-g(3) * (t909 * mrSges(3,1) - t906 * mrSges(3,2)) + t964 * (mrSges(3,1) * t906 + mrSges(3,2) * t909)) * t818 * t916;
t866 = t902 + t899;
t919 = t816 * t923;
t918 = t817 * t923;
t917 = t818 * t923;
t873 = -qJ(3,1) + t883;
t872 = qJ(3,1) + t883;
t871 = -qJ(3,2) + t882;
t870 = qJ(3,2) + t882;
t869 = -qJ(3,3) + t881;
t868 = qJ(3,3) + t881;
t863 = -0.2e1 * qJ(3,1) + t867;
t860 = t915 + t867;
t859 = -0.2e1 * qJ(3,2) + t866;
t856 = t914 + t866;
t855 = -0.2e1 * qJ(3,3) + t865;
t852 = t913 + t865;
t851 = cos(t867);
t850 = cos(t866);
t849 = cos(t865);
t848 = sin(t867);
t847 = sin(t866);
t846 = sin(t865);
t1 = [-t848 * t938 - (t932 * t965 + (cos(t873) + cos(t872)) * t931 + t935 * t958 + (cos(t863) + cos(t860) + 0.2e1 * t851) * pkin(3)) * t917 + (t848 * t957 + cos(t883) * t961 + t851 * t960 - t932 * pkin(3)) * t920 - t847 * t939 - (t933 * t965 + (cos(t871) + cos(t870)) * t931 + t936 * t958 + (cos(t859) + cos(t856) + 0.2e1 * t850) * pkin(3)) * t918 + (t847 * t957 + cos(t882) * t961 + t850 * t960 - t933 * pkin(3)) * t921 - t846 * t940 - (t934 * t965 + (cos(t869) + cos(t868)) * t931 + t937 * t958 + (cos(t855) + cos(t852) + 0.2e1 * t849) * pkin(3)) * t919 + (t846 * t957 + cos(t881) * t961 + t849 * t960 - t934 * pkin(3)) * t922 - g(1) * m(4); t851 * t938 - (t935 * t965 + (sin(t873) + sin(t872)) * t931 + t932 * t957 + (sin(t863) + sin(t860) + 0.2e1 * t848) * pkin(3)) * t917 + (t851 * t958 + sin(t883) * t961 + t848 * t960 - t935 * pkin(3)) * t920 + t850 * t939 - (t936 * t965 + (sin(t871) + sin(t870)) * t931 + t933 * t957 + (sin(t859) + sin(t856) + 0.2e1 * t847) * pkin(3)) * t918 + (t850 * t958 + sin(t882) * t961 + t847 * t960 - t936 * pkin(3)) * t921 + t849 * t940 - (t937 * t965 + (sin(t869) + sin(t868)) * t931 + t934 * t957 + (sin(t855) + sin(t852) + 0.2e1 * t846) * pkin(3)) * t919 + (t849 * t958 + sin(t881) * t961 + t846 * t960 - t937 * pkin(3)) * t922 - g(2) * m(4); (-m(4) - (3 * t912)) * g(3);];
taugX  = t1;
