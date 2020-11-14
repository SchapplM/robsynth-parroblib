% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G1A0
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:01
% EndTime: 2020-08-06 16:50:01
% DurationCPUTime: 0.86s
% Computational Cost: add. (531->138), mult. (1163->263), div. (36->7), fcn. (1074->22), ass. (0->110)
t872 = cos(pkin(3));
t926 = g(3) * t872;
t883 = cos(qJ(3,3));
t878 = sin(qJ(2,3));
t884 = cos(qJ(2,3));
t931 = pkin(2) * t883;
t846 = -pkin(5) * t884 + t878 * t931;
t870 = sin(pkin(3));
t877 = sin(qJ(3,3));
t934 = pkin(2) * t877;
t908 = 0.1e1 / (t846 * t870 + t872 * t934);
t937 = 0.1e1 / t883 * t908;
t885 = cos(qJ(3,2));
t880 = sin(qJ(2,2));
t886 = cos(qJ(2,2));
t930 = pkin(2) * t885;
t847 = -pkin(5) * t886 + t880 * t930;
t879 = sin(qJ(3,2));
t933 = pkin(2) * t879;
t907 = 0.1e1 / (t847 * t870 + t872 * t933);
t936 = 0.1e1 / t885 * t907;
t887 = cos(qJ(3,1));
t882 = sin(qJ(2,1));
t888 = cos(qJ(2,1));
t929 = pkin(2) * t887;
t848 = -pkin(5) * t888 + t882 * t929;
t881 = sin(qJ(3,1));
t932 = pkin(2) * t881;
t906 = 0.1e1 / (t848 * t870 + t872 * t932);
t935 = 0.1e1 / t887 * t906;
t865 = m(1) + m(2) + m(3);
t928 = g(3) * t865;
t927 = g(3) * t870;
t925 = mrSges(3,2) * t870;
t924 = t870 * t883;
t923 = t870 * t885;
t922 = t870 * t887;
t921 = t872 * t878;
t920 = t872 * t880;
t919 = t872 * t882;
t918 = t872 * t884;
t917 = t872 * t886;
t916 = t872 * t888;
t915 = mrSges(3,2) * t926;
t873 = legFrame(3,3);
t859 = sin(t873);
t862 = cos(t873);
t840 = -g(1) * t859 + g(2) * t862;
t843 = g(1) * t862 + g(2) * t859;
t869 = sin(pkin(6));
t871 = cos(pkin(6));
t819 = t840 * t871 - t843 * t869;
t902 = t840 * t869 + t843 * t871;
t896 = t902 * t884;
t905 = -t819 * t872 - t927;
t893 = -t905 * t878 + t896;
t914 = (((t819 * t870 - t926) * mrSges(3,1) + t893 * mrSges(3,2)) * t883 + (t893 * mrSges(3,1) - t819 * t925 + t915) * t877) * t937;
t874 = legFrame(2,3);
t860 = sin(t874);
t863 = cos(t874);
t841 = -g(1) * t860 + g(2) * t863;
t844 = g(1) * t863 + g(2) * t860;
t820 = t841 * t871 - t844 * t869;
t901 = t841 * t869 + t844 * t871;
t895 = t901 * t886;
t904 = -t820 * t872 - t927;
t892 = -t904 * t880 + t895;
t913 = (((t820 * t870 - t926) * mrSges(3,1) + t892 * mrSges(3,2)) * t885 + (t892 * mrSges(3,1) - t820 * t925 + t915) * t879) * t936;
t875 = legFrame(1,3);
t861 = sin(t875);
t864 = cos(t875);
t842 = -g(1) * t861 + g(2) * t864;
t845 = g(1) * t864 + g(2) * t861;
t821 = t842 * t871 - t845 * t869;
t900 = t842 * t869 + t845 * t871;
t894 = t900 * t888;
t903 = -t821 * t872 - t927;
t891 = -t903 * t882 + t894;
t912 = (((t821 * t870 - t926) * mrSges(3,1) + t891 * mrSges(3,2)) * t887 + (t891 * mrSges(3,1) - t821 * t925 + t915) * t881) * t935;
t876 = mrSges(2,2) - mrSges(3,3);
t855 = t876 * t927;
t911 = (t855 * t878 + (t819 * t921 + t896) * t876 + (t902 * t878 + t905 * t884) * (mrSges(3,1) * t883 - mrSges(3,2) * t877 + mrSges(2,1))) * t937;
t910 = (t855 * t880 + (t820 * t920 + t895) * t876 + (t901 * t880 + t904 * t886) * (mrSges(3,1) * t885 - mrSges(3,2) * t879 + mrSges(2,1))) * t936;
t909 = (t855 * t882 + (t821 * t919 + t894) * t876 + (t900 * t882 + t903 * t888) * (mrSges(3,1) * t887 - mrSges(3,2) * t881 + mrSges(2,1))) * t935;
t899 = -t846 * t872 + t870 * t934;
t898 = -t847 * t872 + t870 * t933;
t897 = -t848 * t872 + t870 * t932;
t890 = 0.1e1 / pkin(2);
t851 = pkin(5) * t882 + t888 * t929;
t850 = pkin(5) * t880 + t886 * t930;
t849 = pkin(5) * t878 + t884 * t931;
t839 = -t869 * t919 + t871 * t888;
t838 = -t869 * t920 + t871 * t886;
t837 = -t869 * t921 + t871 * t884;
t836 = t869 * t888 + t871 * t919;
t835 = t869 * t886 + t871 * t920;
t834 = t869 * t884 + t871 * t921;
t833 = t861 * t871 + t864 * t869;
t832 = t860 * t871 + t863 * t869;
t831 = t859 * t871 + t862 * t869;
t830 = -t861 * t869 + t864 * t871;
t829 = -t860 * t869 + t863 * t871;
t828 = -t859 * t869 + t862 * t871;
t818 = t851 * t869 - t897 * t871;
t817 = t850 * t869 - t898 * t871;
t816 = t849 * t869 - t899 * t871;
t815 = t851 * t871 + t897 * t869;
t814 = t850 * t871 + t898 * t869;
t813 = t849 * t871 + t899 * t869;
t1 = [((-t836 * t864 - t839 * t861) * t881 - t830 * t922) * t909 + ((-t835 * t863 - t838 * t860) * t879 - t829 * t923) * t910 + ((-t834 * t862 - t837 * t859) * t877 - t828 * t924) * t911 - g(1) * m(4) + ((-(t830 * t916 - t833 * t882) * t929 - pkin(5) * (t830 * t919 + t833 * t888)) * t912 + (-(t829 * t917 - t832 * t880) * t930 - pkin(5) * (t829 * t920 + t832 * t886)) * t913 + (-(t828 * t918 - t831 * t878) * t931 - pkin(5) * (t828 * t921 + t831 * t884)) * t914) * t890 + (-(t815 * t864 - t818 * t861) * t906 - (t814 * t863 - t817 * t860) * t907 - (t813 * t862 - t816 * t859) * t908) * t928; ((-t836 * t861 + t839 * t864) * t881 - t833 * t922) * t909 + ((-t835 * t860 + t838 * t863) * t879 - t832 * t923) * t910 + ((-t834 * t859 + t837 * t862) * t877 - t831 * t924) * t911 - g(2) * m(4) + ((-(t830 * t882 + t833 * t916) * t929 - (-t830 * t888 + t833 * t919) * pkin(5)) * t912 + (-(t829 * t880 + t832 * t917) * t930 - (-t829 * t886 + t832 * t920) * pkin(5)) * t913 + (-(t828 * t878 + t831 * t918) * t931 - (-t828 * t884 + t831 * t921) * pkin(5)) * t914) * t890 + (-(t815 * t861 + t818 * t864) * t906 - (t814 * t860 + t817 * t863) * t907 - (t813 * t859 + t816 * t862) * t908) * t928; (-m(4) - 0.3e1 * t865) * g(3);];
taugX  = t1;
