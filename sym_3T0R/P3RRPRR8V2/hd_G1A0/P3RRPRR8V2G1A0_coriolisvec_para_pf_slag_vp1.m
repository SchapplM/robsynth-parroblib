% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:10:58
% EndTime: 2022-11-07 13:11:03
% DurationCPUTime: 5.26s
% Computational Cost: add. (26283->432), mult. (28719->675), div. (5061->9), fcn. (20265->62), ass. (0->318)
t826 = rSges(2,3) + pkin(5);
t1004 = 2 * pkin(1);
t961 = qJ(3,1) + pkin(5);
t790 = -pkin(6) - t961;
t779 = 0.1e1 / t790;
t800 = qJ(2,1) + pkin(7);
t758 = sin(t800);
t818 = sin(qJ(2,1));
t970 = pkin(2) * t818;
t723 = pkin(3) * t758 + t970;
t841 = 0.2e1 * qJ(2,1);
t799 = t841 + pkin(7);
t757 = sin(t799);
t995 = pkin(2) * pkin(3);
t928 = 0.2e1 * t995;
t804 = sin(t841);
t849 = pkin(2) ^ 2;
t930 = t849 * t804;
t744 = 0.2e1 * t800;
t736 = sin(t744);
t847 = pkin(3) ^ 2;
t933 = t847 * t736;
t862 = t757 * t928 + t930 + t933;
t692 = t1004 * t723 + t862;
t824 = cos(qJ(2,1));
t783 = t824 * pkin(2);
t769 = cos(t800);
t966 = pkin(3) * t769;
t720 = t783 + t966;
t715 = 0.1e1 / t720;
t831 = xDP(3);
t949 = t715 * t831;
t897 = t692 * t949;
t748 = t783 + pkin(1);
t819 = sin(qJ(1,1));
t825 = cos(qJ(1,1));
t699 = t748 * t819 + t790 * t825;
t702 = t748 * t825 - t790 * t819;
t813 = legFrame(1,3);
t773 = sin(t813);
t776 = cos(t813);
t708 = -t773 * t819 + t776 * t825;
t686 = -t699 * t773 + t702 * t776 + t708 * t966;
t833 = xDP(1);
t955 = t686 * t833;
t707 = t773 * t825 + t776 * t819;
t685 = t699 * t776 + t702 * t773 + t707 * t966;
t832 = xDP(2);
t956 = t685 * t832;
t665 = (t955 + t956 + t897 / 0.2e1) * t779;
t963 = pkin(5) + qJ(3,2);
t789 = -pkin(6) - t963;
t778 = 0.1e1 / t789;
t798 = qJ(2,2) + pkin(7);
t756 = sin(t798);
t816 = sin(qJ(2,2));
t971 = pkin(2) * t816;
t722 = pkin(3) * t756 + t971;
t839 = 0.2e1 * qJ(2,2);
t797 = t839 + pkin(7);
t755 = sin(t797);
t803 = sin(t839);
t931 = t849 * t803;
t743 = 0.2e1 * t798;
t735 = sin(t743);
t934 = t847 * t735;
t863 = t755 * t928 + t931 + t934;
t691 = t1004 * t722 + t863;
t822 = cos(qJ(2,2));
t782 = t822 * pkin(2);
t767 = cos(t798);
t967 = pkin(3) * t767;
t719 = t782 + t967;
t713 = 0.1e1 / t719;
t951 = t713 * t831;
t898 = t691 * t951;
t747 = t782 + pkin(1);
t817 = sin(qJ(1,2));
t823 = cos(qJ(1,2));
t698 = t747 * t817 + t789 * t823;
t701 = t747 * t823 - t789 * t817;
t812 = legFrame(2,3);
t772 = sin(t812);
t775 = cos(t812);
t706 = -t772 * t817 + t775 * t823;
t684 = -t698 * t772 + t701 * t775 + t706 * t967;
t957 = t684 * t833;
t705 = t772 * t823 + t775 * t817;
t683 = t698 * t775 + t701 * t772 + t705 * t967;
t958 = t683 * t832;
t664 = (t957 + t958 + t898 / 0.2e1) * t778;
t962 = pkin(5) + qJ(3,3);
t788 = -pkin(6) - t962;
t777 = 0.1e1 / t788;
t796 = qJ(2,3) + pkin(7);
t754 = sin(t796);
t814 = sin(qJ(2,3));
t972 = pkin(2) * t814;
t721 = pkin(3) * t754 + t972;
t837 = 0.2e1 * qJ(2,3);
t795 = t837 + pkin(7);
t753 = sin(t795);
t802 = sin(t837);
t932 = t849 * t802;
t742 = 0.2e1 * t796;
t734 = sin(t742);
t935 = t847 * t734;
t864 = t753 * t928 + t932 + t935;
t690 = t1004 * t721 + t864;
t820 = cos(qJ(2,3));
t781 = t820 * pkin(2);
t764 = cos(t796);
t968 = pkin(3) * t764;
t718 = t781 + t968;
t711 = 0.1e1 / t718;
t953 = t711 * t831;
t899 = t690 * t953;
t746 = t781 + pkin(1);
t815 = sin(qJ(1,3));
t821 = cos(qJ(1,3));
t697 = t746 * t815 + t788 * t821;
t700 = t746 * t821 - t788 * t815;
t811 = legFrame(3,3);
t771 = sin(t811);
t774 = cos(t811);
t704 = -t771 * t815 + t774 * t821;
t682 = -t697 * t771 + t700 * t774 + t704 * t968;
t959 = t682 * t833;
t703 = t771 * t821 + t774 * t815;
t681 = t697 * t774 + t700 * t771 + t703 * t968;
t960 = t681 * t832;
t663 = (t959 + t960 + t899 / 0.2e1) * t777;
t1001 = 0.1e1 / t718 ^ 2;
t678 = (t703 * t832 + t704 * t833 + t721 * t953) * t777;
t675 = pkin(1) * t678;
t654 = -t675 - (-t959 / 0.2e1 - t960 / 0.2e1 - t899 / 0.4e1) * t777;
t850 = pkin(1) ^ 2;
t669 = (t788 ^ 2 + t850) * t678;
t687 = t788 * t953;
t835 = 0.2e1 * pkin(7);
t759 = cos(t835 + qJ(2,3));
t765 = cos(qJ(2,3) - pkin(7));
t808 = t831 ^ 2;
t836 = 0.3e1 * qJ(2,3);
t846 = pkin(3) * t847;
t848 = pkin(2) * t849;
t737 = cos(t742);
t784 = t847 + t849;
t805 = cos(t837);
t869 = t737 * t847 + t805 * t849 + t784;
t883 = 0.3e1 / 0.4e1 * t849;
t884 = 0.3e1 / 0.4e1 * t847;
t965 = pkin(3) * t849;
t885 = -0.2e1 * t846 - 0.4e1 * t965;
t954 = t711 * t777;
t888 = -t954 / 0.2e1;
t910 = -0.2e1 * t965;
t969 = pkin(2) * t847;
t911 = -0.2e1 * t969;
t919 = -0.6e1 * t849 - 0.3e1 * t847;
t763 = cos(t795);
t810 = cos(pkin(7));
t923 = t763 + t810;
t983 = -0.3e1 / 0.4e1 * t849;
t984 = -t831 / 0.2e1;
t991 = -t678 / 0.4e1;
t997 = -0.2e1 * t848 - 0.4e1 * t969;
t964 = t810 * pkin(2);
t904 = pkin(3) * t964;
t998 = -0.4e1 * pkin(1) * (t904 + t849 / 0.2e1 + t847 / 0.2e1);
t645 = (t759 * t911 + t764 * t885 + t765 * t910 + t820 * t997 + t998) * t808 * t1001 * t888 + (t864 * t1001 * t788 * t984 + ((-t846 * cos(0.3e1 * t796) - t848 * cos(t836)) * t991 - (t935 / 0.2e1 + t932 / 0.2e1) * t687 + (-(-cos(t836 + pkin(7)) - t765) * t678 * t883 + (-t1004 * t663 + t919 * t991 + t669) * t764) * pkin(3) + ((-t678 * t983 + t669) * t820 - (-cos(t835 + t836) - t759 - 0.2e1 * t820) * t678 * t884 + (-0.2e1 * t654 * t923 - t687 * t753) * pkin(3) - (pkin(3) * t923 + t1004 * t820) * t663) * pkin(2) + (-t663 / 0.2e1 - t654) * t869) * t711) * t777 * t678;
t712 = t711 * t1001;
t890 = -0.2e1 * t904;
t902 = t678 * t954;
t929 = -0.2e1 * t995;
t945 = (0.2e1 * t904 + t784) * t808;
t648 = t712 * t777 * t945 + (-(t763 * t929 - t869 + t890) * t678 / 0.2e1 + (-0.2e1 * t663 + t675) * t718) * t902;
t785 = rSges(3,3) + t962;
t981 = m(3) * t785;
t728 = rSges(3,2) * t981 - Icges(3,6);
t731 = rSges(3,1) * t981 - Icges(3,5);
t809 = sin(pkin(7));
t994 = m(2) * rSges(2,1);
t873 = -t826 * t994 + Icges(2,5);
t993 = m(2) * rSges(2,2);
t889 = t826 * t993 - Icges(2,6);
t907 = pkin(2) * t981;
t672 = -(t728 * t809 - t731 * t810 + t873 - t907) * t814 + (t728 * t810 + t731 * t809 + t889) * t820;
t694 = -rSges(3,1) * t764 + rSges(3,2) * t754 - t746;
t982 = m(2) * t826;
t740 = rSges(2,2) * t982 - Icges(2,6);
t992 = m(3) * rSges(3,2);
t741 = t809 * pkin(2) * t992;
t834 = 2 * t850;
t843 = (rSges(2,2) ^ 2);
t845 = (rSges(2,1) ^ 2);
t920 = (t843 + t845);
t865 = -((2 * pkin(5) ^ 2 + t834 + (4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3) + t920) * m(2)) / 0.2e1 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t741 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t842 = rSges(3,2) ^ 2;
t844 = rSges(3,1) ^ 2;
t866 = -t834 / 0.2e1 - t842 / 0.2e1 - t844 / 0.2e1 - t849 / 0.2e1;
t996 = pkin(2) * m(3);
t745 = t994 + t996;
t872 = t745 * t814 + t820 * t993;
t895 = t712 * t721 * t808;
t903 = t993 * t1004;
t908 = 0.2e1 * rSges(3,2);
t909 = 0.2e1 * rSges(3,1);
t912 = -0.2e1 * pkin(1) * t745;
t1003 = 0.2e1 * m(3);
t915 = t678 * t1003;
t917 = rSges(3,2) * t1004;
t927 = rSges(2,1) * t982 - Icges(2,5);
t750 = rSges(2,1) * t993 - Icges(2,4);
t938 = t750 * t805;
t749 = -rSges(3,1) * t992 + Icges(3,4);
t941 = t749 * t737;
t726 = m(3) * (-t842 + t844) - Icges(3,1) + Icges(3,2);
t944 = t726 * t734;
t717 = t849 * m(3) + ((-t843 + t845) * m(2)) + Icges(2,2) - Icges(2,1);
t948 = t717 * t802;
t975 = pkin(2) * t753;
t978 = pkin(1) * t764;
t985 = t750 / 0.2e1;
t986 = -t749 / 0.2e1;
t987 = -t726 / 0.2e1;
t988 = -t717 / 0.2e1;
t926 = (t737 * t987 + t805 * t988 + t814 * t903 + t820 * t912 + t865) * t648 - t672 * t895 + 0.2e1 * (t734 * t986 + t802 * t985) * t648 + (-t694 * t645 + (rSges(3,2) * t975 + t754 * t917 - t785 ^ 2 + t866 + (-pkin(2) * t923 - 0.2e1 * t978) * rSges(3,1)) * t648) * m(3) + t663 * t785 * t915 + ((-t731 * t764 + t728 * t754 - (t907 + t927) * t820 + t740 * t814) * t953 - (-t948 - t944 + 0.2e1 * t941 - 0.2e1 * t938 - t872 * t1004 + ((-pkin(2) * t763 - t978) * t908 + (-pkin(1) * t754 - t975) * t909) * m(3)) * t678) * t953;
t1000 = 0.1e1 / t719 ^ 2;
t679 = (t705 * t832 + t706 * t833 + t722 * t951) * t778;
t676 = pkin(1) * t679;
t655 = -t676 - (-t957 / 0.2e1 - t958 / 0.2e1 - t898 / 0.4e1) * t778;
t670 = (t789 ^ 2 + t850) * t679;
t688 = t789 * t951;
t760 = cos(t835 + qJ(2,2));
t762 = cos(-pkin(7) + qJ(2,2));
t838 = 0.3e1 * qJ(2,2);
t738 = cos(t743);
t806 = cos(t839);
t868 = t847 * t738 + t849 * t806 + t784;
t952 = t713 * t778;
t887 = -t952 / 0.2e1;
t766 = cos(t797);
t922 = t766 + t810;
t990 = -t679 / 0.4e1;
t646 = (t760 * t911 + t762 * t910 + t767 * t885 + t822 * t997 + t998) * t808 * t1000 * t887 + (t863 * t1000 * t789 * t984 + ((-t846 * cos(0.3e1 * t798) - t848 * cos(t838)) * t990 - (t934 / 0.2e1 + t931 / 0.2e1) * t688 + (-(-cos(pkin(7) + t838) - t762) * t679 * t883 + (-t1004 * t664 + t919 * t990 + t670) * t767) * pkin(3) + ((-t679 * t983 + t670) * t822 - (-cos(t835 + t838) - t760 - 0.2e1 * t822) * t679 * t884 + (-0.2e1 * t655 * t922 - t688 * t755) * pkin(3) - (pkin(3) * t922 + t1004 * t822) * t664) * pkin(2) + (-t655 - t664 / 0.2e1) * t868) * t713) * t778 * t679;
t714 = t713 * t1000;
t901 = t679 * t952;
t649 = t714 * t778 * t945 + (-(t766 * t929 - t868 + t890) * t679 / 0.2e1 + (-0.2e1 * t664 + t676) * t719) * t901;
t786 = rSges(3,3) + t963;
t980 = m(3) * t786;
t729 = rSges(3,2) * t980 - Icges(3,6);
t732 = rSges(3,1) * t980 - Icges(3,5);
t906 = pkin(2) * t980;
t673 = -(t729 * t809 - t732 * t810 + t873 - t906) * t816 + (t729 * t810 + t732 * t809 + t889) * t822;
t695 = -rSges(3,1) * t767 + rSges(3,2) * t756 - t747;
t871 = t745 * t816 + t822 * t993;
t893 = t714 * t722 * t808;
t914 = t679 * t1003;
t937 = t750 * t806;
t940 = t749 * t738;
t943 = t726 * t735;
t947 = t717 * t803;
t974 = pkin(2) * t755;
t977 = pkin(1) * t767;
t925 = (t738 * t987 + t806 * t988 + t816 * t903 + t822 * t912 + t865) * t649 - t673 * t893 + 0.2e1 * (t735 * t986 + t803 * t985) * t649 + (-t695 * t646 + (rSges(3,2) * t974 + t756 * t917 - t786 ^ 2 + t866 + (-pkin(2) * t922 - 0.2e1 * t977) * rSges(3,1)) * t649) * m(3) + t664 * t786 * t914 + ((-t732 * t767 + t729 * t756 - (t906 + t927) * t822 + t740 * t816) * t951 - (-t947 - t943 + 0.2e1 * t940 - 0.2e1 * t937 - t871 * t1004 + ((-pkin(2) * t766 - t977) * t908 + (-pkin(1) * t756 - t974) * t909) * m(3)) * t679) * t951;
t680 = (t707 * t832 + t708 * t833 + t723 * t949) * t779;
t677 = pkin(1) * t680;
t656 = -t677 - (-t955 / 0.2e1 - t956 / 0.2e1 - t897 / 0.4e1) * t779;
t671 = (t790 ^ 2 + t850) * t680;
t689 = t790 * t949;
t761 = cos(t835 + qJ(2,1));
t770 = cos(qJ(2,1) - pkin(7));
t840 = 0.3e1 * qJ(2,1);
t739 = cos(t744);
t807 = cos(t841);
t867 = t847 * t739 + t849 * t807 + t784;
t950 = t715 * t779;
t886 = -t950 / 0.2e1;
t768 = cos(t799);
t921 = t768 + t810;
t989 = -t680 / 0.4e1;
t999 = 0.1e1 / t720 ^ 2;
t647 = (t761 * t911 + t769 * t885 + t770 * t910 + t824 * t997 + t998) * t808 * t999 * t886 + (t862 * t999 * t790 * t984 + ((-t846 * cos(0.3e1 * t800) - t848 * cos(t840)) * t989 - (t933 / 0.2e1 + t930 / 0.2e1) * t689 + (-(-cos(t840 + pkin(7)) - t770) * t680 * t883 + (-t1004 * t665 + t919 * t989 + t671) * t769) * pkin(3) + ((-t680 * t983 + t671) * t824 - (-cos(t840 + t835) - t761 - 0.2e1 * t824) * t680 * t884 + (-0.2e1 * t656 * t921 - t689 * t757) * pkin(3) - (pkin(3) * t921 + t1004 * t824) * t665) * pkin(2) + (-t656 - t665 / 0.2e1) * t867) * t715) * t779 * t680;
t716 = t715 * t999;
t900 = t680 * t950;
t650 = t716 * t779 * t945 + (-(t768 * t929 - t867 + t890) * t680 / 0.2e1 + (-0.2e1 * t665 + t677) * t720) * t900;
t787 = rSges(3,3) + t961;
t979 = m(3) * t787;
t730 = rSges(3,2) * t979 - Icges(3,6);
t733 = rSges(3,1) * t979 - Icges(3,5);
t905 = pkin(2) * t979;
t674 = -(t730 * t809 - t733 * t810 + t873 - t905) * t818 + (t730 * t810 + t733 * t809 + t889) * t824;
t696 = -rSges(3,1) * t769 + rSges(3,2) * t758 - t748;
t870 = t745 * t818 + t824 * t993;
t891 = t716 * t723 * t808;
t913 = t680 * t1003;
t936 = t750 * t807;
t939 = t749 * t739;
t942 = t726 * t736;
t946 = t717 * t804;
t973 = pkin(2) * t757;
t976 = pkin(1) * t769;
t924 = (t739 * t987 + t807 * t988 + t818 * t903 + t824 * t912 + t865) * t650 - t674 * t891 + 0.2e1 * (t736 * t986 + t804 * t985) * t650 + (-t696 * t647 + (rSges(3,2) * t973 + t758 * t917 - t787 ^ 2 + t866 + (-pkin(2) * t921 - 0.2e1 * t976) * rSges(3,1)) * t650) * m(3) + t665 * t787 * t913 + ((-t733 * t769 + t730 * t758 - (t905 + t927) * t824 + t740 * t818) * t949 - (-t946 - t942 + 0.2e1 * t939 - 0.2e1 * t936 - t870 * t1004 + ((-pkin(2) * t768 - t976) * t908 + (-pkin(1) * t758 - t973) * t909) * m(3)) * t680) * t949;
t642 = (-t648 * t694 - t645) * m(3);
t876 = rSges(3,1) * t754 + rSges(3,2) * t764;
t666 = t678 * t785 / 0.2e1 + (t876 + t972) * t953;
t882 = -t666 * t915 + t642;
t643 = (-t649 * t695 - t646) * m(3);
t875 = rSges(3,1) * t756 + rSges(3,2) * t767;
t667 = t679 * t786 / 0.2e1 + (t875 + t971) * t951;
t881 = -t667 * t914 + t643;
t644 = (-t650 * t696 - t647) * m(3);
t874 = rSges(3,1) * t758 + rSges(3,2) * t769;
t668 = t680 * t787 / 0.2e1 + (t874 + t970) * t949;
t880 = -t668 * t913 + t644;
t693 = 0.2e1 * t741 - (t920 * m(2)) - Icges(2,3) - Icges(3,3) + (-0.2e1 * rSges(3,1) * t964 - t842 - t844 - t849) * m(3);
t1 = [-(t686 * t880 + t708 * t924) * t779 - (t684 * t881 + t706 * t925) * t778 - (t682 * t882 + t704 * t926) * t777; -(t685 * t880 + t707 * t924) * t779 - (t683 * t881 + t705 * t925) * t778 - (t681 * t882 + t703 * t926) * t777; -t926 * t721 * t954 - t925 * t722 * t952 - t924 * t723 * t950 + (t674 * t650 - t693 * t891 - ((0.2e1 * t665 * t970 + t874 * (-t677 - (-t897 - 0.2e1 * t955 - 0.2e1 * t956) * t779)) * m(3) - (t942 / 0.2e1 - t939 + t946 / 0.2e1 + t936 + t870 * pkin(1) + (rSges(3,1) * t757 + rSges(3,2) * t768) * t996) * t680) * t680) * t715 + (t673 * t649 - t693 * t893 - ((0.2e1 * t664 * t971 + t875 * (-t676 - (-t898 - 0.2e1 * t957 - 0.2e1 * t958) * t778)) * m(3) - (t943 / 0.2e1 - t940 + t947 / 0.2e1 + t937 + t871 * pkin(1) + (rSges(3,1) * t755 + rSges(3,2) * t766) * t996) * t679) * t679) * t713 + (t672 * t648 - t693 * t895 - ((0.2e1 * t663 * t972 + t876 * (-t675 - (-t899 - 0.2e1 * t959 - 0.2e1 * t960) * t777)) * m(3) - (t944 / 0.2e1 - t941 + t948 / 0.2e1 + t938 + t872 * pkin(1) + (rSges(3,1) * t753 + rSges(3,2) * t763) * t996) * t678) * t678) * t711 + (m(3) * t668 * t900 + t644 * t886) * t692 + (m(3) * t667 * t901 + t643 * t887) * t691 + (m(3) * t666 * t902 + t642 * t888) * t690;];
taucX  = t1;
