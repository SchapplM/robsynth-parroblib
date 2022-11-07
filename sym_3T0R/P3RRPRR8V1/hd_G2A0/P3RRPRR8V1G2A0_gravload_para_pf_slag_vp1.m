% Calculate Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:17
% EndTime: 2022-11-04 17:04:18
% DurationCPUTime: 0.48s
% Computational Cost: add. (564->124), mult. (1095->201), div. (30->9), fcn. (618->23), ass. (0->89)
t910 = sin(pkin(5));
t973 = cos(pkin(5));
t882 = -m(2) * rSges(2,1) + (-rSges(3,1) * t973 + rSges(3,2) * t910 - pkin(1)) * m(3);
t883 = m(2) * rSges(2,2) + (rSges(3,1) * t910 + t973 * rSges(3,2)) * m(3);
t924 = sin(qJ(2,1));
t930 = cos(qJ(2,1));
t981 = -t930 * t882 - t883 * t924;
t922 = sin(qJ(2,2));
t928 = cos(qJ(2,2));
t980 = -t928 * t882 - t883 * t922;
t920 = sin(qJ(2,3));
t926 = cos(qJ(2,3));
t979 = -t926 * t882 - t883 * t920;
t978 = m(1) * rSges(1,1);
t977 = m(3) * (rSges(3,3) + qJ(3,3));
t976 = m(3) * (rSges(3,3) + qJ(3,2));
t975 = m(3) * (rSges(3,3) + qJ(3,1));
t974 = pkin(2) * t910;
t917 = legFrame(3,2);
t900 = sin(t917);
t903 = cos(t917);
t887 = t903 * g(1) - t900 * g(2);
t914 = pkin(4) + qJ(3,3);
t907 = 0.1e1 / t914;
t921 = sin(qJ(1,3));
t927 = cos(qJ(1,3));
t972 = (-g(3) * t921 + t887 * t927) * t907;
t918 = legFrame(2,2);
t901 = sin(t918);
t904 = cos(t918);
t888 = t904 * g(1) - t901 * g(2);
t915 = pkin(4) + qJ(3,2);
t908 = 0.1e1 / t915;
t923 = sin(qJ(1,2));
t929 = cos(qJ(1,2));
t971 = (-g(3) * t923 + t888 * t929) * t908;
t919 = legFrame(1,2);
t902 = sin(t919);
t905 = cos(t919);
t889 = t905 * g(1) - t902 * g(2);
t916 = pkin(4) + qJ(3,1);
t909 = 0.1e1 / t916;
t925 = sin(qJ(1,1));
t931 = cos(qJ(1,1));
t970 = (-g(3) * t925 + t889 * t931) * t909;
t884 = t900 * g(1) + t903 * g(2);
t939 = g(3) * t927 + t887 * t921;
t942 = t926 * pkin(1) + pkin(2) * cos(qJ(2,3) + pkin(5));
t966 = 0.1e1 / t942 * ((-t939 * t882 + t884 * t883) * t920 + t926 * (t884 * t882 + t939 * t883));
t885 = t901 * g(1) + t904 * g(2);
t938 = g(3) * t929 + t888 * t923;
t941 = t928 * pkin(1) + pkin(2) * cos(qJ(2,2) + pkin(5));
t965 = 0.1e1 / t941 * ((-t938 * t882 + t885 * t883) * t922 + t928 * (t885 * t882 + t938 * t883));
t886 = t902 * g(1) + t905 * g(2);
t937 = g(3) * t931 + t889 * t925;
t940 = t930 * pkin(1) + pkin(2) * cos(qJ(2,1) + pkin(5));
t964 = 0.1e1 / t940 * ((-t937 * t882 + t886 * t883) * t924 + t930 * (t886 * t882 + t937 * t883));
t894 = t973 * pkin(2) + pkin(1);
t963 = t894 * t903;
t962 = t894 * t904;
t961 = t894 * t905;
t960 = t900 * t894;
t959 = t901 * t894;
t958 = t902 * t894;
t895 = m(1) * rSges(1,2) - m(2) * rSges(2,3);
t893 = t895 * g(3);
t906 = g(3) * t978;
t957 = t907 * ((-g(3) * t977 + t893) * t927 + t921 * (t979 * g(3) + t906) + ((-t978 - t979) * t927 + t921 * (t895 - t977)) * t887);
t956 = t908 * ((-g(3) * t976 + t893) * t929 + t923 * (t980 * g(3) + t906) + ((-t978 - t980) * t929 + t923 * (t895 - t976)) * t888);
t955 = t909 * ((-g(3) * t975 + t893) * t931 + t925 * (t981 * g(3) + t906) + ((-t978 - t981) * t931 + t925 * (t895 - t975)) * t889);
t951 = t903 * t974;
t950 = t904 * t974;
t949 = t905 * t974;
t948 = t900 * t974;
t947 = t901 * t974;
t946 = t902 * t974;
t936 = t894 * t926 - t920 * t974;
t945 = 0.1e1 / t936 * t957;
t935 = t894 * t928 - t922 * t974;
t944 = 0.1e1 / t935 * t956;
t934 = t894 * t930 - t924 * t974;
t943 = 0.1e1 / t934 * t955;
t881 = t924 * t894 + t930 * t974;
t880 = t922 * t894 + t928 * t974;
t879 = t920 * t894 + t926 * t974;
t872 = -t931 * t916 + t934 * t925;
t871 = -t929 * t915 + t935 * t923;
t870 = -t927 * t914 + t936 * t921;
t1 = [((t925 * t961 + t946) * t930 + (-t925 * t949 + t958) * t924) * t943 + t902 * t964 + ((t923 * t962 + t947) * t928 + (-t923 * t950 + t959) * t922) * t944 + t901 * t965 + ((t921 * t963 + t948) * t926 + (-t921 * t951 + t960) * t920) * t945 + t900 * t966 - m(4) * g(1) + ((t872 * t905 + t902 * t881) * t970 + (t871 * t904 + t901 * t880) * t971 + (t870 * t903 + t900 * t879) * t972) * m(3); ((-t925 * t958 + t949) * t930 + (t925 * t946 + t961) * t924) * t943 + t905 * t964 + ((-t923 * t959 + t950) * t928 + (t923 * t947 + t962) * t922) * t944 + t904 * t965 + ((-t921 * t960 + t951) * t926 + (t921 * t948 + t963) * t920) * t945 + t903 * t966 - m(4) * g(2) + ((-t872 * t902 + t905 * t881) * t970 + (-t871 * t901 + t904 * t880) * t971 + (-t870 * t900 + t903 * t879) * t972) * m(3); t927 * t957 + t929 * t956 + t931 * t955 - m(4) * g(3) + ((t925 * t916 + t940 * t931) * t970 + (t923 * t915 + t941 * t929) * t971 + (t921 * t914 + t942 * t927) * t972) * m(3);];
taugX  = t1;
