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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
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

function taugX = P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:56
% EndTime: 2020-08-06 18:31:57
% DurationCPUTime: 0.94s
% Computational Cost: add. (722->195), mult. (703->224), div. (24->7), fcn. (462->82), ass. (0->111)
t971 = 2 * pkin(2);
t970 = -2 * pkin(1);
t969 = -2 * pkin(2);
t928 = (-pkin(6) - pkin(5));
t967 = -2 * t928;
t966 = 2 * t928;
t965 = rSges(3,3) + pkin(5);
t964 = m(3) / pkin(3);
t916 = sin(qJ(1,3));
t963 = t916 * pkin(1);
t918 = sin(qJ(1,2));
t962 = t918 * pkin(1);
t920 = sin(qJ(1,1));
t961 = t920 * pkin(1);
t960 = cos(pkin(7)) * pkin(1) + pkin(2);
t912 = legFrame(3,3);
t900 = sin(t912);
t903 = cos(t912);
t841 = -t900 * g(1) + t903 * g(2);
t844 = t903 * g(1) + t900 * g(2);
t909 = qJ(1,3) + pkin(7);
t890 = sin(t909);
t893 = cos(t909);
t922 = cos(qJ(1,3));
t906 = t922 * pkin(1);
t921 = cos(qJ(3,3));
t915 = sin(qJ(3,3));
t941 = t921 * rSges(3,1) - t915 * rSges(3,2);
t935 = pkin(2) + t941;
t959 = 0.1e1 / (t921 * pkin(3) + t960) * (-m(1) * (t844 * (-t916 * rSges(1,1) - t922 * rSges(1,2)) + t841 * (t922 * rSges(1,1) - t916 * rSges(1,2))) - m(2) * (t844 * (-t890 * rSges(2,1) - t893 * rSges(2,2) - t963) + t841 * (t893 * rSges(2,1) - t890 * rSges(2,2) + t906)) - m(3) * (-t844 * t963 + t841 * t906 + (t841 * t935 + t844 * t965) * t893 + (t841 * t965 - t844 * t935) * t890));
t913 = legFrame(2,3);
t901 = sin(t913);
t904 = cos(t913);
t842 = -t901 * g(1) + t904 * g(2);
t845 = t904 * g(1) + t901 * g(2);
t910 = qJ(1,2) + pkin(7);
t891 = sin(t910);
t894 = cos(t910);
t924 = cos(qJ(1,2));
t907 = t924 * pkin(1);
t923 = cos(qJ(3,2));
t917 = sin(qJ(3,2));
t939 = t923 * rSges(3,1) - t917 * rSges(3,2);
t934 = pkin(2) + t939;
t958 = 0.1e1 / (t923 * pkin(3) + t960) * (-m(1) * (t845 * (-t918 * rSges(1,1) - t924 * rSges(1,2)) + t842 * (t924 * rSges(1,1) - t918 * rSges(1,2))) - m(2) * (t845 * (-t891 * rSges(2,1) - t894 * rSges(2,2) - t962) + t842 * (t894 * rSges(2,1) - t891 * rSges(2,2) + t907)) - m(3) * (-t845 * t962 + t842 * t907 + (t842 * t934 + t845 * t965) * t894 + (t842 * t965 - t845 * t934) * t891));
t914 = legFrame(1,3);
t902 = sin(t914);
t905 = cos(t914);
t843 = -t902 * g(1) + t905 * g(2);
t846 = t905 * g(1) + t902 * g(2);
t911 = qJ(1,1) + pkin(7);
t892 = sin(t911);
t895 = cos(t911);
t926 = cos(qJ(1,1));
t908 = t926 * pkin(1);
t925 = cos(qJ(3,1));
t919 = sin(qJ(3,1));
t937 = t925 * rSges(3,1) - t919 * rSges(3,2);
t933 = pkin(2) + t937;
t957 = 0.1e1 / (t925 * pkin(3) + t960) * (-m(1) * (t846 * (-t920 * rSges(1,1) - t926 * rSges(1,2)) + t843 * (t926 * rSges(1,1) - t920 * rSges(1,2))) - m(2) * (t846 * (-t892 * rSges(2,1) - t895 * rSges(2,2) - t961) + t843 * (t895 * rSges(2,1) - t892 * rSges(2,2) + t908)) - m(3) * (-t846 * t961 + t843 * t908 + (t843 * t933 + t846 * t965) * t895 + (t843 * t965 - t846 * t933) * t892));
t956 = (-m(2) - m(3)) * g(3) / 0.2e1;
t898 = qJ(1,1) + t914;
t897 = qJ(1,2) + t913;
t896 = qJ(1,3) + t912;
t880 = t912 + t909;
t869 = qJ(3,3) + t880;
t870 = -qJ(3,3) + t880;
t955 = sin(t869) + sin(t870);
t881 = t913 + t910;
t873 = qJ(3,2) + t881;
t874 = -qJ(3,2) + t881;
t954 = sin(t873) + sin(t874);
t882 = t914 + t911;
t877 = qJ(3,1) + t882;
t878 = -qJ(3,1) + t882;
t953 = sin(t877) + sin(t878);
t952 = cos(t869) + cos(t870);
t951 = cos(t873) + cos(t874);
t950 = cos(t877) + cos(t878);
t949 = 2 * pkin(1);
t929 = 0.2e1 * qJ(3,3);
t838 = 0.1e1 / (pkin(3) * sin(t929) + t915 * t971 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1));
t947 = t838 * t956;
t930 = 0.2e1 * qJ(3,2);
t839 = 0.1e1 / (pkin(3) * sin(t930) + t917 * t971 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1));
t946 = t839 * t956;
t931 = 0.2e1 * qJ(3,1);
t840 = 0.1e1 / (pkin(3) * sin(t931) + t919 * t971 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1));
t945 = t840 * t956;
t944 = (g(3) * t941 + (t841 * t890 + t844 * t893) * (-rSges(3,1) * t915 - rSges(3,2) * t921)) * t838 * t964;
t943 = (g(3) * t939 + (t842 * t891 + t845 * t894) * (-rSges(3,1) * t917 - rSges(3,2) * t923)) * t839 * t964;
t942 = (g(3) * t937 + (t843 * t892 + t846 * t895) * (-rSges(3,1) * t919 - rSges(3,2) * t925)) * t840 * t964;
t888 = -qJ(3,1) + t898;
t887 = qJ(3,1) + t898;
t886 = -qJ(3,2) + t897;
t885 = qJ(3,2) + t897;
t884 = -qJ(3,3) + t896;
t883 = qJ(3,3) + t896;
t879 = -0.2e1 * qJ(3,1) + t882;
t876 = t931 + t882;
t875 = -0.2e1 * qJ(3,2) + t881;
t872 = t930 + t881;
t871 = -0.2e1 * qJ(3,3) + t880;
t868 = t929 + t880;
t867 = cos(t882);
t866 = cos(t881);
t865 = cos(t880);
t864 = sin(t882);
t863 = sin(t881);
t862 = sin(t880);
t1 = [-t864 * t957 + (t950 * t971 + (cos(t888) + cos(t887)) * t949 + t953 * t967 + (cos(t879) + cos(t876) + 0.2e1 * t867) * pkin(3)) * t945 - (t864 * t966 + cos(t898) * t970 + t867 * t969 - t950 * pkin(3)) * t942 - t863 * t958 + (t951 * t971 + (cos(t886) + cos(t885)) * t949 + t954 * t967 + (cos(t875) + cos(t872) + 0.2e1 * t866) * pkin(3)) * t946 - (t863 * t966 + cos(t897) * t970 + t866 * t969 - t951 * pkin(3)) * t943 - t862 * t959 + (t952 * t971 + (cos(t884) + cos(t883)) * t949 + t955 * t967 + (cos(t871) + cos(t868) + 0.2e1 * t865) * pkin(3)) * t947 - (t862 * t966 + cos(t896) * t970 + t865 * t969 - t952 * pkin(3)) * t944 - m(4) * g(1); t867 * t957 + (t953 * t971 + (sin(t888) + sin(t887)) * t949 + t950 * t966 + (sin(t879) + sin(t876) + 0.2e1 * t864) * pkin(3)) * t945 - (t867 * t967 + sin(t898) * t970 + t864 * t969 - t953 * pkin(3)) * t942 + t866 * t958 + (t954 * t971 + (sin(t886) + sin(t885)) * t949 + t951 * t966 + (sin(t875) + sin(t872) + 0.2e1 * t863) * pkin(3)) * t946 - (t866 * t967 + sin(t897) * t970 + t863 * t969 - t954 * pkin(3)) * t943 + t865 * t959 + (t955 * t971 + (sin(t884) + sin(t883)) * t949 + t952 * t966 + (sin(t871) + sin(t868) + 0.2e1 * t862) * pkin(3)) * t947 - (t865 * t967 + sin(t896) * t970 + t862 * t969 - t955 * pkin(3)) * t944 - m(4) * g(2); (-0.3e1 * m(2) - 0.3e1 * m(3) - m(4)) * g(3);];
taugX  = t1;
