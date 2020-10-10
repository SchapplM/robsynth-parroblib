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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
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

function taugX = P4PRRR1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:12:39
% EndTime: 2020-03-02 20:12:39
% DurationCPUTime: 0.83s
% Computational Cost: add. (910->119), mult. (819->202), div. (80->6), fcn. (714->26), ass. (0->101)
t919 = pkin(7) + qJ(2,4);
t898 = qJ(3,4) + t919;
t888 = sin(t898);
t889 = cos(t898);
t896 = sin(t919);
t897 = cos(t919);
t951 = 0.1e1 / (t888 * t897 - t896 * t889);
t920 = pkin(7) + qJ(2,3);
t905 = qJ(3,3) + t920;
t890 = sin(t905);
t893 = cos(t905);
t899 = sin(t920);
t902 = cos(t920);
t950 = 0.1e1 / (t890 * t902 - t899 * t893);
t921 = pkin(7) + qJ(2,2);
t906 = qJ(3,2) + t921;
t891 = sin(t906);
t894 = cos(t906);
t900 = sin(t921);
t903 = cos(t921);
t949 = 0.1e1 / (t891 * t903 - t900 * t894);
t922 = pkin(7) + qJ(2,1);
t907 = qJ(3,1) + t922;
t892 = sin(t907);
t895 = cos(t907);
t901 = sin(t922);
t904 = cos(t922);
t948 = 0.1e1 / (t892 * t904 - t901 * t895);
t923 = legFrame(4,3);
t908 = sin(t923);
t912 = cos(t923);
t880 = -t908 * g(1) + t912 * g(2);
t884 = t912 * g(1) + t908 * g(2);
t832 = -t889 * (mrSges(3,1) * t880 - mrSges(3,2) * t884) + (mrSges(3,1) * t884 + mrSges(3,2) * t880) * t888;
t916 = m(3) * pkin(2) + mrSges(2,1);
t947 = ((mrSges(2,2) * t884 - t880 * t916) * t897 + (t880 * mrSges(2,2) + t916 * t884) * t896 + t832) * t951;
t924 = legFrame(3,3);
t909 = sin(t924);
t913 = cos(t924);
t881 = -t909 * g(1) + t913 * g(2);
t885 = t913 * g(1) + t909 * g(2);
t833 = -t893 * (mrSges(3,1) * t881 - mrSges(3,2) * t885) + (mrSges(3,1) * t885 + mrSges(3,2) * t881) * t890;
t946 = ((mrSges(2,2) * t885 - t881 * t916) * t902 + (t881 * mrSges(2,2) + t916 * t885) * t899 + t833) * t950;
t925 = legFrame(2,3);
t910 = sin(t925);
t914 = cos(t925);
t882 = -t910 * g(1) + t914 * g(2);
t886 = t914 * g(1) + t910 * g(2);
t834 = -t894 * (mrSges(3,1) * t882 - mrSges(3,2) * t886) + (mrSges(3,1) * t886 + mrSges(3,2) * t882) * t891;
t945 = ((mrSges(2,2) * t886 - t882 * t916) * t903 + (t882 * mrSges(2,2) + t916 * t886) * t900 + t834) * t949;
t926 = legFrame(1,3);
t911 = sin(t926);
t915 = cos(t926);
t883 = -t911 * g(1) + t915 * g(2);
t887 = t915 * g(1) + t911 * g(2);
t835 = -t895 * (mrSges(3,1) * t883 - mrSges(3,2) * t887) + (mrSges(3,1) * t887 + mrSges(3,2) * t883) * t892;
t944 = ((mrSges(2,2) * t887 - t883 * t916) * t904 + (t883 * mrSges(2,2) + t916 * t887) * t901 + t835) * t948;
t943 = t832 * t951;
t942 = t833 * t950;
t941 = t834 * t949;
t940 = t835 * t948;
t864 = t912 * t888 + t908 * t889;
t865 = -t888 * t908 + t912 * t889;
t866 = t913 * t890 + t909 * t893;
t867 = -t890 * t909 + t913 * t893;
t868 = t914 * t891 + t910 * t894;
t869 = -t891 * t910 + t914 * t894;
t870 = t915 * t892 + t911 * t895;
t871 = -t892 * t911 + t915 * t895;
t939 = 0.1e1 / pkin(2);
t938 = 0.1e1 / pkin(3);
t937 = koppelP(1,1);
t936 = koppelP(2,1);
t935 = koppelP(3,1);
t934 = koppelP(4,1);
t933 = koppelP(1,2);
t932 = koppelP(2,2);
t931 = koppelP(3,2);
t930 = koppelP(4,2);
t929 = mrSges(4,1);
t928 = mrSges(4,2);
t927 = xP(4);
t918 = cos(t927);
t917 = sin(t927);
t879 = -t917 * t933 + t918 * t937;
t878 = -t917 * t932 + t918 * t936;
t877 = -t917 * t931 + t918 * t935;
t876 = -t917 * t930 + t918 * t934;
t875 = -t917 * t937 - t918 * t933;
t874 = -t917 * t936 - t918 * t932;
t873 = -t917 * t935 - t918 * t931;
t872 = -t917 * t934 - t918 * t930;
t843 = -pkin(2) * (t901 * t911 - t915 * t904) + t871 * pkin(3);
t842 = -pkin(2) * (t900 * t910 - t914 * t903) + t869 * pkin(3);
t841 = -pkin(2) * (t899 * t909 - t913 * t902) + t867 * pkin(3);
t840 = pkin(2) * (t915 * t901 + t911 * t904) + t870 * pkin(3);
t839 = pkin(2) * (t914 * t900 + t910 * t903) + t868 * pkin(3);
t838 = pkin(2) * (t913 * t899 + t909 * t902) + t866 * pkin(3);
t837 = -pkin(2) * (t896 * t908 - t912 * t897) + t865 * pkin(3);
t836 = pkin(2) * (t912 * t896 + t908 * t897) + t864 * pkin(3);
t1 = [-g(1) * m(4) + (t865 * t947 + t867 * t946 + t869 * t945 + t871 * t944 + (-t837 * t943 - t841 * t942 - t842 * t941 - t843 * t940) * t938) * t939; -g(2) * m(4) + (t864 * t947 + t866 * t946 + t868 * t945 + t870 * t944 + (-t836 * t943 - t838 * t942 - t839 * t941 - t840 * t940) * t938) * t939; -(4 * g(3) * (m(1) + m(2) + m(3))) - g(3) * m(4); -(-g(1) * t929 - g(2) * t928) * t917 + t918 * (g(1) * t928 - g(2) * t929) + ((t870 * t879 + t871 * t875) * t944 + (t868 * t878 + t869 * t874) * t945 + (t866 * t877 + t867 * t873) * t946 + (t864 * t876 + t865 * t872) * t947 + (-(t840 * t879 + t843 * t875) * t940 - (t839 * t878 + t842 * t874) * t941 - (t838 * t877 + t841 * t873) * t942 - (t836 * t876 + t837 * t872) * t943) * t938) * t939;];
taugX  = t1;
