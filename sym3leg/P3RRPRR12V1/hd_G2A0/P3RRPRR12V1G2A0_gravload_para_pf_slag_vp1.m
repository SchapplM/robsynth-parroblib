% Calculate Gravitation load for parallel robot
% P3RRPRR12V1G2A0
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
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:36
% EndTime: 2020-08-06 19:05:37
% DurationCPUTime: 0.66s
% Computational Cost: add. (591->122), mult. (954->213), div. (36->6), fcn. (540->18), ass. (0->98)
t915 = m(2) * rSges(2,2);
t882 = (-qJ(3,1) - rSges(3,3)) * m(3) + t915;
t886 = (pkin(1) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t907 = sin(qJ(2,1));
t913 = cos(qJ(2,1));
t962 = -t882 * t907 + t886 * t913;
t881 = (-qJ(3,2) - rSges(3,3)) * m(3) + t915;
t905 = sin(qJ(2,2));
t911 = cos(qJ(2,2));
t961 = -t881 * t905 + t886 * t911;
t880 = (-qJ(3,3) - rSges(3,3)) * m(3) + t915;
t903 = sin(qJ(2,3));
t909 = cos(qJ(2,3));
t960 = -t880 * t903 + t886 * t909;
t959 = m(1) * rSges(1,1);
t910 = cos(qJ(1,3));
t958 = t910 * pkin(4);
t912 = cos(qJ(1,2));
t957 = t912 * pkin(4);
t914 = cos(qJ(1,1));
t956 = t914 * pkin(4);
t900 = legFrame(3,2);
t887 = sin(t900);
t890 = cos(t900);
t860 = t887 * g(1) + t890 * g(2);
t917 = 0.1e1 / qJ(3,3);
t863 = t890 * g(1) - t887 * g(2);
t904 = sin(qJ(1,3));
t922 = g(3) * t910 + t863 * t904;
t955 = ((-t886 * t860 + t922 * t880) * t909 + (t860 * t880 + t922 * t886) * t903) * t917;
t901 = legFrame(2,2);
t888 = sin(t901);
t891 = cos(t901);
t861 = t888 * g(1) + t891 * g(2);
t918 = 0.1e1 / qJ(3,2);
t864 = t891 * g(1) - t888 * g(2);
t906 = sin(qJ(1,2));
t921 = g(3) * t912 + t864 * t906;
t954 = ((-t886 * t861 + t921 * t881) * t911 + (t861 * t881 + t921 * t886) * t905) * t918;
t902 = legFrame(1,2);
t889 = sin(t902);
t892 = cos(t902);
t862 = t889 * g(1) + t892 * g(2);
t919 = 0.1e1 / qJ(3,1);
t865 = t892 * g(1) - t889 * g(2);
t908 = sin(qJ(1,1));
t920 = g(3) * t914 + t865 * t908;
t953 = ((-t886 * t862 + t920 * t882) * t913 + (t862 * t882 + t920 * t886) * t907) * t919;
t952 = (-t860 * t909 + t922 * t903) * t917;
t951 = (-t861 * t911 + t921 * t905) * t918;
t950 = (-t862 * t913 + t920 * t907) * t919;
t943 = t887 * qJ(3,3);
t942 = t888 * qJ(3,2);
t941 = t889 * qJ(3,1);
t940 = t890 * qJ(3,3);
t939 = t891 * qJ(3,2);
t938 = t892 * qJ(3,1);
t937 = t903 * qJ(3,3);
t916 = pkin(1) + pkin(2);
t936 = t904 * t916;
t935 = t905 * qJ(3,2);
t934 = t906 * t916;
t933 = t907 * qJ(3,1);
t932 = t908 * t916;
t879 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - rSges(3,2) * m(3);
t878 = t879 * g(3);
t893 = g(3) * t959;
t848 = t878 * t910 + (t960 * g(3) + t893) * t904 + ((-t959 - t960) * t910 + t879 * t904) * t863;
t931 = t910 * t848;
t849 = t878 * t912 + (t961 * g(3) + t893) * t906 + ((-t959 - t961) * t912 + t879 * t906) * t864;
t930 = t912 * t849;
t850 = t878 * t914 + (t962 * g(3) + t893) * t908 + ((-t959 - t962) * t914 + t879 * t908) * t865;
t929 = t914 * t850;
t928 = t916 * t903;
t927 = t916 * t905;
t926 = t916 * t907;
t875 = t916 * t909 + t937;
t925 = -t904 * pkin(4) + t875 * t910;
t876 = t916 * t911 + t935;
t924 = -t906 * pkin(4) + t876 * t912;
t877 = t916 * t913 + t933;
t923 = -t908 * pkin(4) + t877 * t914;
t896 = t913 ^ 2;
t895 = t911 ^ 2;
t894 = t909 ^ 2;
t874 = -t913 * qJ(3,1) + t926;
t873 = -t911 * qJ(3,2) + t927;
t872 = -t909 * qJ(3,3) + t928;
t871 = 0.1e1 / t877;
t870 = 0.1e1 / t876;
t869 = 0.1e1 / t875;
t868 = t908 * t933 + t956;
t867 = t906 * t935 + t957;
t866 = t904 * t937 + t958;
t859 = t877 * t908 + t956;
t858 = t876 * t906 + t957;
t857 = t875 * t904 + t958;
t1 = [-m(4) * g(1) + (t892 * t929 + ((t892 * t932 - t941) * t896 + (t868 * t892 + t889 * t926) * t913 + t941) * t953) * t871 + (t891 * t930 + ((t891 * t934 - t942) * t895 + (t867 * t891 + t888 * t927) * t911 + t942) * t954) * t870 + (t890 * t931 + ((t890 * t936 - t943) * t894 + (t866 * t890 + t887 * t928) * t909 + t943) * t955) * t869 + (-(t859 * t892 + t889 * t874) * t950 - (t858 * t891 + t888 * t873) * t951 - (t857 * t890 + t887 * t872) * t952) * m(3); -m(4) * g(2) + (-t889 * t929 + ((-t889 * t932 - t938) * t896 + (-t889 * t868 + t892 * t926) * t913 + t938) * t953) * t871 + (-t888 * t930 + ((-t888 * t934 - t939) * t895 + (-t888 * t867 + t891 * t927) * t911 + t939) * t954) * t870 + (-t887 * t931 + ((-t887 * t936 - t940) * t894 + (-t887 * t866 + t890 * t928) * t909 + t940) * t955) * t869 + (-(-t859 * t889 + t874 * t892) * t950 - (-t858 * t888 + t873 * t891) * t951 - (-t857 * t887 + t872 * t890) * t952) * m(3); -m(4) * g(3) + (t913 * t923 * t953 - t908 * t850) * t871 + (t911 * t924 * t954 - t906 * t849) * t870 + (t909 * t925 * t955 - t904 * t848) * t869 + (-t923 * t950 - t924 * t951 - t925 * t952) * m(3);];
taugX  = t1;
