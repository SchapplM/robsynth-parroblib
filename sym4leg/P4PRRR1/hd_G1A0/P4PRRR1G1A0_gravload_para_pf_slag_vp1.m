% Calculate Gravitation load for parallel robot
% P4PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRR1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:10:39
% EndTime: 2020-03-02 20:10:40
% DurationCPUTime: 0.95s
% Computational Cost: add. (934->126), mult. (952->222), div. (80->6), fcn. (762->26), ass. (0->100)
t928 = pkin(7) + qJ(2,4);
t908 = qJ(3,4) + t928;
t898 = sin(t908);
t899 = cos(t908);
t906 = sin(t928);
t907 = cos(t928);
t961 = 0.1e1 / (t898 * t907 - t906 * t899);
t929 = pkin(7) + qJ(2,3);
t915 = qJ(3,3) + t929;
t900 = sin(t915);
t903 = cos(t915);
t909 = sin(t929);
t912 = cos(t929);
t960 = 0.1e1 / (t900 * t912 - t909 * t903);
t930 = pkin(7) + qJ(2,2);
t916 = qJ(3,2) + t930;
t901 = sin(t916);
t904 = cos(t916);
t910 = sin(t930);
t913 = cos(t930);
t959 = 0.1e1 / (t901 * t913 - t910 * t904);
t931 = pkin(7) + qJ(2,1);
t917 = qJ(3,1) + t931;
t902 = sin(t917);
t905 = cos(t917);
t911 = sin(t931);
t914 = cos(t931);
t958 = 0.1e1 / (t902 * t914 - t911 * t905);
t957 = m(3) / pkin(3);
t932 = legFrame(4,3);
t918 = sin(t932);
t922 = cos(t932);
t890 = -t918 * g(1) + t922 * g(2);
t894 = t922 * g(1) + t918 * g(2);
t842 = -(rSges(3,1) * t890 - rSges(3,2) * t894) * t899 + (rSges(3,1) * t894 + rSges(3,2) * t890) * t898;
t956 = (((-rSges(2,1) * t890 + rSges(2,2) * t894) * t907 + (rSges(2,1) * t894 + rSges(2,2) * t890) * t906) * m(2) + ((-t890 * t907 + t894 * t906) * pkin(2) + t842) * m(3)) * t961;
t933 = legFrame(3,3);
t919 = sin(t933);
t923 = cos(t933);
t891 = -t919 * g(1) + t923 * g(2);
t895 = t923 * g(1) + t919 * g(2);
t843 = -(rSges(3,1) * t891 - rSges(3,2) * t895) * t903 + (rSges(3,1) * t895 + rSges(3,2) * t891) * t900;
t955 = (((-rSges(2,1) * t891 + rSges(2,2) * t895) * t912 + (rSges(2,1) * t895 + rSges(2,2) * t891) * t909) * m(2) + ((-t891 * t912 + t895 * t909) * pkin(2) + t843) * m(3)) * t960;
t934 = legFrame(2,3);
t920 = sin(t934);
t924 = cos(t934);
t892 = -t920 * g(1) + t924 * g(2);
t896 = t924 * g(1) + t920 * g(2);
t844 = -(rSges(3,1) * t892 - rSges(3,2) * t896) * t904 + (rSges(3,1) * t896 + rSges(3,2) * t892) * t901;
t954 = (((-rSges(2,1) * t892 + rSges(2,2) * t896) * t913 + (rSges(2,1) * t896 + rSges(2,2) * t892) * t910) * m(2) + ((-t892 * t913 + t896 * t910) * pkin(2) + t844) * m(3)) * t959;
t935 = legFrame(1,3);
t921 = sin(t935);
t925 = cos(t935);
t893 = -t921 * g(1) + t925 * g(2);
t897 = t925 * g(1) + t921 * g(2);
t845 = -(rSges(3,1) * t893 - rSges(3,2) * t897) * t905 + (rSges(3,1) * t897 + rSges(3,2) * t893) * t902;
t953 = (((-rSges(2,1) * t893 + rSges(2,2) * t897) * t914 + (rSges(2,1) * t897 + rSges(2,2) * t893) * t911) * m(2) + ((-t893 * t914 + t897 * t911) * pkin(2) + t845) * m(3)) * t958;
t952 = t842 * t961;
t951 = t843 * t960;
t950 = t844 * t959;
t949 = t845 * t958;
t874 = t922 * t898 + t918 * t899;
t875 = -t898 * t918 + t922 * t899;
t876 = t923 * t900 + t919 * t903;
t877 = -t900 * t919 + t923 * t903;
t878 = t924 * t901 + t920 * t904;
t879 = -t901 * t920 + t924 * t904;
t880 = t925 * t902 + t921 * t905;
t881 = -t902 * t921 + t925 * t905;
t948 = 0.1e1 / pkin(2);
t946 = koppelP(1,1);
t945 = koppelP(2,1);
t944 = koppelP(3,1);
t943 = koppelP(4,1);
t942 = koppelP(1,2);
t941 = koppelP(2,2);
t940 = koppelP(3,2);
t939 = koppelP(4,2);
t938 = rSges(4,1);
t937 = rSges(4,2);
t936 = xP(4);
t927 = cos(t936);
t926 = sin(t936);
t889 = -t926 * t942 + t927 * t946;
t888 = -t926 * t941 + t927 * t945;
t887 = -t926 * t940 + t927 * t944;
t886 = -t926 * t939 + t927 * t943;
t885 = -t926 * t946 - t927 * t942;
t884 = -t926 * t945 - t927 * t941;
t883 = -t926 * t944 - t927 * t940;
t882 = -t926 * t943 - t927 * t939;
t853 = -pkin(2) * (t911 * t921 - t925 * t914) + t881 * pkin(3);
t852 = -pkin(2) * (t910 * t920 - t924 * t913) + t879 * pkin(3);
t851 = -pkin(2) * (t909 * t919 - t923 * t912) + t877 * pkin(3);
t850 = pkin(2) * (t925 * t911 + t921 * t914) + t880 * pkin(3);
t849 = pkin(2) * (t924 * t910 + t920 * t913) + t878 * pkin(3);
t848 = pkin(2) * (t923 * t909 + t919 * t912) + t876 * pkin(3);
t847 = -pkin(2) * (t906 * t918 - t922 * t907) + t875 * pkin(3);
t846 = pkin(2) * (t922 * t906 + t918 * t907) + t874 * pkin(3);
t1 = [-m(4) * g(1) + (t875 * t956 + t877 * t955 + t879 * t954 + t881 * t953 + (-t847 * t952 - t851 * t951 - t852 * t950 - t853 * t949) * t957) * t948; -m(4) * g(2) + (t874 * t956 + t876 * t955 + t878 * t954 + t880 * t953 + (-t846 * t952 - t848 * t951 - t849 * t950 - t850 * t949) * t957) * t948; -0.4e1 * g(3) * (m(1) + m(2) + m(3)) - m(4) * g(3); ((g(1) * t937 - g(2) * t938) * t927 + (g(1) * t938 + g(2) * t937) * t926) * m(4) + ((t880 * t889 + t881 * t885) * t953 + (t878 * t888 + t879 * t884) * t954 + (t876 * t887 + t877 * t883) * t955 + (t874 * t886 + t875 * t882) * t956 + (-(t850 * t889 + t853 * t885) * t949 - (t849 * t888 + t852 * t884) * t950 - (t848 * t887 + t851 * t883) * t951 - (t846 * t886 + t847 * t882) * t952) * t957) * t948;];
taugX  = t1;
