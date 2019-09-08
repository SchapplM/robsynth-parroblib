% Calculate Gravitation load for parallel robot
% P3RPP1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPP1G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:06
% EndTime: 2019-05-03 14:52:06
% DurationCPUTime: 0.56s
% Computational Cost: add. (765->139), mult. (1014->216), div. (36->3), fcn. (452->14), ass. (0->131)
t916 = 2 * pkin(1);
t915 = pkin(1) ^ 2 + 1;
t870 = pkin(1) + qJ(3,1);
t869 = pkin(1) + qJ(3,2);
t868 = pkin(1) + qJ(3,3);
t899 = (pkin(1) - rSges(2,2)) * m(2) + m(1) * rSges(1,1);
t823 = (rSges(3,3) + t868) * m(3) + t899;
t865 = legFrame(3,3);
t850 = sin(t865);
t853 = cos(t865);
t826 = -t850 * g(1) + t853 * g(2);
t829 = t853 * g(1) + t850 * g(2);
t877 = m(1) * rSges(1,2);
t832 = (-qJ(2,3) - rSges(3,2)) * m(3) + (-qJ(2,3) - rSges(2,3)) * m(2) + t877;
t871 = sin(qJ(1,3));
t874 = cos(qJ(1,3));
t778 = (-t826 * t823 + t832 * t829) * t874 + (t823 * t829 + t826 * t832) * t871;
t844 = t915 + (t916 + qJ(3,3)) * qJ(3,3);
t884 = qJ(2,3) ^ 2;
t841 = 1 / (t884 + t844);
t914 = t778 * t841;
t824 = (rSges(3,3) + t869) * m(3) + t899;
t866 = legFrame(2,3);
t851 = sin(t866);
t854 = cos(t866);
t827 = -t851 * g(1) + t854 * g(2);
t830 = t854 * g(1) + t851 * g(2);
t833 = (-qJ(2,2) - rSges(3,2)) * m(3) + (-qJ(2,2) - rSges(2,3)) * m(2) + t877;
t872 = sin(qJ(1,2));
t875 = cos(qJ(1,2));
t779 = (-t827 * t824 + t833 * t830) * t875 + (t824 * t830 + t827 * t833) * t872;
t845 = t915 + (t916 + qJ(3,2)) * qJ(3,2);
t886 = qJ(2,2) ^ 2;
t842 = 1 / (t886 + t845);
t913 = t779 * t842;
t825 = (rSges(3,3) + t870) * m(3) + t899;
t867 = legFrame(1,3);
t852 = sin(t867);
t855 = cos(t867);
t828 = -t852 * g(1) + t855 * g(2);
t831 = t855 * g(1) + t852 * g(2);
t834 = (-qJ(2,1) - rSges(3,2)) * m(3) + (-qJ(2,1) - rSges(2,3)) * m(2) + t877;
t873 = sin(qJ(1,1));
t876 = cos(qJ(1,1));
t780 = (-t828 * t825 + t834 * t831) * t876 + (t825 * t831 + t828 * t834) * t873;
t846 = t915 + (t916 + qJ(3,1)) * qJ(3,1);
t888 = qJ(2,1) ^ 2;
t843 = 1 / (t888 + t846);
t912 = t780 * t843;
t799 = t826 * t874 - t829 * t871;
t911 = t799 * t841;
t800 = t827 * t875 - t830 * t872;
t910 = t800 * t842;
t801 = t828 * t876 - t831 * t873;
t909 = t801 * t843;
t802 = t826 * t871 + t829 * t874;
t908 = t802 * t841;
t803 = t827 * t872 + t830 * t875;
t907 = t803 * t842;
t804 = t828 * t873 + t831 * t876;
t906 = t804 * t843;
t905 = t871 * t868;
t904 = t872 * t869;
t903 = t873 * t870;
t902 = t874 * qJ(2,3);
t901 = t875 * qJ(2,2);
t900 = t876 * qJ(2,1);
t898 = t870 * t900;
t897 = t869 * t901;
t896 = t868 * t902;
t894 = koppelP(1,1);
t893 = koppelP(2,1);
t892 = koppelP(3,1);
t891 = koppelP(1,2);
t890 = koppelP(2,2);
t889 = koppelP(3,2);
t882 = rSges(4,1);
t881 = rSges(4,2);
t880 = xP(3);
t879 = m(2) + m(3);
t864 = 1 + t888;
t863 = 1 + t886;
t862 = 1 + t884;
t858 = cos(t880);
t857 = sin(t880);
t849 = qJ(2,1) * t903;
t848 = qJ(2,2) * t904;
t847 = qJ(2,3) * t905;
t840 = t873 * qJ(2,1) + t870 * t876;
t839 = t872 * qJ(2,2) + t869 * t875;
t838 = t871 * qJ(2,3) + t868 * t874;
t837 = t900 - t903;
t836 = t901 - t904;
t835 = t902 - t905;
t822 = -t857 * t891 + t858 * t894;
t821 = -t857 * t890 + t858 * t893;
t820 = -t857 * t889 + t858 * t892;
t819 = -t857 * t894 - t858 * t891;
t818 = -t857 * t893 - t858 * t890;
t817 = -t857 * t892 - t858 * t889;
t816 = -t864 * t876 + t849;
t815 = -t863 * t875 + t848;
t814 = -t862 * t874 + t847;
t813 = t873 * t864 + t898;
t812 = t872 * t863 + t897;
t811 = t871 * t862 + t896;
t810 = t846 * t876 + t849;
t809 = t845 * t875 + t848;
t808 = t844 * t874 + t847;
t807 = t873 * t846 - t898;
t806 = t872 * t845 - t897;
t805 = t871 * t844 - t896;
t798 = t852 * t837 + t840 * t855;
t797 = t851 * t836 + t839 * t854;
t796 = t850 * t835 + t838 * t853;
t795 = t837 * t855 - t852 * t840;
t794 = t836 * t854 - t851 * t839;
t793 = t835 * t853 - t850 * t838;
t792 = t852 * t813 + t816 * t855;
t791 = t851 * t812 + t815 * t854;
t790 = t850 * t811 + t814 * t853;
t789 = t813 * t855 - t852 * t816;
t788 = t812 * t854 - t851 * t815;
t787 = t811 * t853 - t850 * t814;
t786 = -t852 * t807 + t810 * t855;
t785 = -t851 * t806 + t809 * t854;
t784 = -t850 * t805 + t808 * t853;
t783 = t807 * t855 + t852 * t810;
t782 = t806 * t854 + t851 * t809;
t781 = t805 * t853 + t850 * t808;
t1 = [t793 * t914 + t794 * t913 + t795 * t912 - m(4) * g(1) + (t787 * t911 + t788 * t910 + t789 * t909) * t879 + (-t784 * t908 - t785 * t907 - t786 * t906) * m(3); t796 * t914 + t797 * t913 + t798 * t912 - m(4) * g(2) + (t790 * t911 + t791 * t910 + t792 * t909) * t879 + (-t781 * t908 - t782 * t907 - t783 * t906) * m(3); m(4) * ((g(1) * t882 + g(2) * t881) * t857 + (g(1) * t881 - g(2) * t882) * t858) + ((t795 * t819 + t798 * t822) * t780 + (t789 * t819 + t792 * t822) * t801 * t879 - (t783 * t822 + t786 * t819) * m(3) * t804) * t843 + ((t794 * t818 + t797 * t821) * t779 + (t788 * t818 + t791 * t821) * t800 * t879 - (t782 * t821 + t785 * t818) * m(3) * t803) * t842 + ((t793 * t817 + t796 * t820) * t778 + (t787 * t817 + t790 * t820) * t799 * t879 - (t781 * t820 + t784 * t817) * m(3) * t802) * t841;];
taugX  = t1;
