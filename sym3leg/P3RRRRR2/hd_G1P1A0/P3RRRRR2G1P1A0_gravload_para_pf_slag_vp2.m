% Calculate Gravitation load for parallel robot
% P3RRRRR2G1P1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:20
% EndTime: 2020-03-09 21:04:20
% DurationCPUTime: 0.45s
% Computational Cost: add. (729->147), mult. (1020->193), div. (57->14), fcn. (663->72), ass. (0->103)
t960 = -2 * pkin(1);
t914 = legFrame(3,3);
t893 = sin(t914);
t896 = cos(t914);
t862 = -t893 * g(1) + t896 * g(2);
t959 = mrSges(3,1) * t862;
t915 = legFrame(2,3);
t894 = sin(t915);
t897 = cos(t915);
t863 = -t894 * g(1) + t897 * g(2);
t958 = mrSges(3,1) * t863;
t916 = legFrame(1,3);
t895 = sin(t916);
t898 = cos(t916);
t864 = -t895 * g(1) + t898 * g(2);
t957 = mrSges(3,1) * t864;
t956 = mrSges(3,2) * t862;
t955 = mrSges(3,2) * t863;
t954 = mrSges(3,2) * t864;
t865 = t896 * g(1) + t893 * g(2);
t953 = mrSges(3,2) * t865;
t866 = t897 * g(1) + t894 * g(2);
t952 = mrSges(3,2) * t866;
t867 = t898 * g(1) + t895 * g(2);
t951 = mrSges(3,2) * t867;
t925 = 1 / pkin(1);
t950 = 0.1e1 / sin(qJ(2,3)) * t925;
t949 = 0.1e1 / sin(qJ(2,2)) * t925;
t948 = 0.1e1 / sin(qJ(2,1)) * t925;
t911 = qJ(1,3) + qJ(2,3);
t890 = cos(t911);
t917 = mrSges(3,3) - mrSges(2,2);
t947 = t917 * t890;
t912 = qJ(2,2) + qJ(1,2);
t891 = cos(t912);
t946 = t917 * t891;
t913 = qJ(1,1) + qJ(2,1);
t892 = cos(t913);
t945 = t917 * t892;
t924 = 1 / pkin(2);
t944 = t924 * t925;
t886 = qJ(1,1) + t916;
t943 = qJ(2,1) - qJ(3,1);
t942 = qJ(2,1) + qJ(3,1);
t941 = qJ(2,2) - qJ(3,2);
t940 = qJ(2,2) + qJ(3,2);
t884 = qJ(1,3) + t914;
t939 = qJ(2,3) - qJ(3,3);
t938 = qJ(2,3) + qJ(3,3);
t856 = mrSges(3,1) * t865;
t859 = mrSges(2,1) * t865;
t883 = (m(2) + m(3)) * pkin(1) + mrSges(1,1);
t887 = sin(t911);
t899 = qJ(1,3) + t938;
t900 = qJ(1,3) + t939;
t850 = (-t953 - t959) * cos(t900) / 0.2e1 + (t856 - t956) * sin(t900) / 0.2e1 + (t953 - t959) * cos(t899) / 0.2e1 + (t856 + t956) * sin(t899) / 0.2e1 - mrSges(2,1) * t862 * t890 + (t862 * mrSges(1,2) + t883 * t865) * sin(qJ(1,3)) - t865 * t947 + (-t862 * t917 + t859) * t887 + (mrSges(1,2) * t865 - t883 * t862) * cos(qJ(1,3));
t937 = t850 * t950;
t857 = mrSges(3,1) * t866;
t860 = mrSges(2,1) * t866;
t888 = sin(t912);
t901 = qJ(1,2) + t940;
t902 = qJ(1,2) + t941;
t851 = (-t952 - t958) * cos(t902) / 0.2e1 + (t857 - t955) * sin(t902) / 0.2e1 + (t952 - t958) * cos(t901) / 0.2e1 + (t857 + t955) * sin(t901) / 0.2e1 - mrSges(2,1) * t863 * t891 + (t863 * mrSges(1,2) + t883 * t866) * sin(qJ(1,2)) - t866 * t946 + (-t863 * t917 + t860) * t888 + (mrSges(1,2) * t866 - t883 * t863) * cos(qJ(1,2));
t936 = t851 * t949;
t858 = mrSges(3,1) * t867;
t861 = mrSges(2,1) * t867;
t889 = sin(t913);
t903 = qJ(1,1) + t942;
t904 = qJ(1,1) + t943;
t852 = (-t951 - t957) * cos(t904) / 0.2e1 + (t858 - t954) * sin(t904) / 0.2e1 + (t951 - t957) * cos(t903) / 0.2e1 + (t858 + t954) * sin(t903) / 0.2e1 - mrSges(2,1) * t864 * t892 + (t864 * mrSges(1,2) + t883 * t867) * sin(qJ(1,1)) - t867 * t945 + (-t864 * t917 + t861) * t889 + (mrSges(1,2) * t867 - t883 * t864) * cos(qJ(1,1));
t935 = t852 * t948;
t918 = sin(qJ(3,3));
t934 = t918 * t950;
t919 = sin(qJ(3,2));
t933 = t919 * t949;
t920 = sin(qJ(3,1));
t932 = t920 * t948;
t882 = qJ(2,1) + t886;
t880 = qJ(2,3) + t884;
t921 = cos(qJ(3,3));
t928 = t921 * mrSges(3,1) - t918 * mrSges(3,2);
t853 = t859 * t887 + (t928 * t887 - t947) * t865 + ((-mrSges(2,1) - t928) * t890 - t917 * t887) * t862;
t931 = t853 / (sin(t938) + sin(t939)) * t944;
t922 = cos(qJ(3,2));
t927 = t922 * mrSges(3,1) - t919 * mrSges(3,2);
t854 = t860 * t888 + (t927 * t888 - t946) * t866 + ((-mrSges(2,1) - t927) * t891 - t917 * t888) * t863;
t930 = t854 / (sin(t940) + sin(t941)) * t944;
t923 = cos(qJ(3,1));
t926 = t923 * mrSges(3,1) - t920 * mrSges(3,2);
t855 = t861 * t889 + (t926 * t889 - t945) * t867 + ((-mrSges(2,1) - t926) * t892 - t917 * t889) * t864;
t929 = t855 / (sin(t942) + sin(t943)) * t944;
t910 = 0.1e1 / t923;
t909 = 0.1e1 / t922;
t908 = 0.1e1 / t921;
t885 = qJ(1,2) + t915;
t881 = t915 + t912;
t879 = -qJ(3,1) + t882;
t878 = qJ(3,1) + t882;
t877 = t915 + t902;
t876 = t915 + t901;
t875 = -qJ(3,3) + t880;
t874 = qJ(3,3) + t880;
t1 = [cos(t882) * t935 + (cos(t886) * t960 + (-cos(t879) - cos(t878)) * pkin(2)) * t929 + cos(t881) * t936 + (cos(t885) * t960 + (-cos(t877) - cos(t876)) * pkin(2)) * t930 + cos(t880) * t937 + (cos(t884) * t960 + (-cos(t875) - cos(t874)) * pkin(2)) * t931 - g(1) * m(4); sin(t882) * t935 + (sin(t886) * t960 + (-sin(t879) - sin(t878)) * pkin(2)) * t929 + sin(t881) * t936 + (sin(t885) * t960 + (-sin(t877) - sin(t876)) * pkin(2)) * t930 + sin(t880) * t937 + (sin(t884) * t960 + (-sin(t875) - sin(t874)) * pkin(2)) * t931 - g(2) * m(4); t908 * t850 * t934 + t909 * t851 * t933 + t910 * t852 * t932 - g(3) * m(4) + (-(cos(qJ(2,1)) * pkin(1) + t923 * pkin(2)) / t923 ^ 2 * t855 * t932 + t910 * (-g(3) * t926 + (t864 * t889 + t867 * t892) * (mrSges(3,1) * t920 + mrSges(3,2) * t923)) - (cos(qJ(2,2)) * pkin(1) + t922 * pkin(2)) / t922 ^ 2 * t854 * t933 + t909 * (-g(3) * t927 + (t863 * t888 + t866 * t891) * (mrSges(3,1) * t919 + mrSges(3,2) * t922)) - (cos(qJ(2,3)) * pkin(1) + t921 * pkin(2)) / t921 ^ 2 * t853 * t934 + t908 * (-g(3) * t928 + (t862 * t887 + t865 * t890) * (mrSges(3,1) * t918 + mrSges(3,2) * t921))) * t924;];
taugX  = t1;
