% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR6V1G1A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:31:53
% EndTime: 2020-08-06 18:31:56
% DurationCPUTime: 2.30s
% Computational Cost: add. (12390->303), mult. (11061->443), div. (1395->13), fcn. (6093->68), ass. (0->204)
t694 = cos(pkin(7)) * pkin(1);
t716 = cos(qJ(3,3));
t805 = t716 * pkin(3) + pkin(2);
t630 = t694 + t805;
t627 = 0.1e1 / t630;
t688 = qJ(1,3) + legFrame(3,3);
t676 = pkin(7) + t688;
t655 = sin(t676);
t658 = cos(t676);
t724 = xDP(2);
t725 = xDP(1);
t609 = (-t655 * t725 + t658 * t724) * t627;
t726 = -pkin(6) - pkin(5);
t806 = pkin(1) * sin(pkin(7));
t819 = -t726 + t806;
t831 = t609 * t819;
t717 = cos(qJ(3,2));
t804 = t717 * pkin(3) + pkin(2);
t631 = t694 + t804;
t628 = 0.1e1 / t631;
t689 = qJ(1,2) + legFrame(2,3);
t677 = pkin(7) + t689;
t656 = sin(t677);
t659 = cos(t677);
t610 = (-t656 * t725 + t659 * t724) * t628;
t830 = t610 * t819;
t718 = cos(qJ(3,1));
t803 = t718 * pkin(3) + pkin(2);
t632 = t694 + t803;
t629 = 0.1e1 / t632;
t690 = qJ(1,1) + legFrame(1,3);
t678 = pkin(7) + t690;
t657 = sin(t678);
t660 = cos(t678);
t611 = (-t657 * t725 + t660 * t724) * t629;
t829 = t611 * t819;
t737 = pkin(1) ^ 2;
t828 = -pkin(2) ^ 2 - t737;
t817 = -0.2e1 * pkin(2);
t827 = 0.2e1 * pkin(2);
t826 = m(3) / 0.2e1;
t670 = qJ(3,1) + t678;
t671 = -qJ(3,1) + t678;
t772 = sin(t670) + sin(t671);
t813 = -0.2e1 * t726;
t818 = -0.2e1 * pkin(1);
t602 = t660 * t813 + sin(t690) * t818 + t657 * t817 - t772 * pkin(3);
t769 = cos(t670) + cos(t671);
t812 = 0.2e1 * t726;
t605 = t657 * t812 + cos(t690) * t818 + t660 * t817 - t769 * pkin(3);
t731 = 0.2e1 * qJ(3,1);
t700 = sin(t731);
t693 = pkin(3) * t700;
t715 = sin(qJ(3,1));
t617 = 0.1e1 / (t715 * t827 + t693 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1));
t735 = 0.1e1 / pkin(3);
t590 = (t602 * t724 + t605 * t725) * t735 * t617;
t587 = t590 ^ 2;
t608 = t611 ^ 2;
t685 = t694 + pkin(2);
t777 = t715 * t685;
t621 = 0.1e1 / (t693 + 0.2e1 * t777);
t703 = cos(t731);
t706 = t718 ^ 2;
t776 = t715 * t718;
t814 = -0.2e1 * t685;
t734 = pkin(3) ^ 2;
t747 = t726 ^ 2 + t806 * t813 - t828;
t763 = 0.2e1 * t694;
t815 = t763 * t817 - t734 - 0.2e1 * t747;
t578 = (t608 * t718 * t815 + (-0.2e1 * t587 * t632 + ((-t632 * t703 + t706 * t814 - t685) * t611 + (t700 + 0.2e1 * t776) * t819 * t590) * t611) * pkin(3)) * t621;
t728 = -m(2) - m(3);
t575 = t728 * t578;
t783 = t590 * t715;
t759 = pkin(3) * t783;
t775 = pkin(3) * t827;
t782 = t590 * t718;
t581 = (-t819 * t759 + (t706 * t734 + t718 * t775 + t803 * t763 + t747) * t611) / (t693 / 0.2e1 + t777) * t611 * t735 + 0.2e1 * (t632 * t782 - t715 * t829) * t621 * t590;
t639 = rSges(3,1) * t718 - rSges(3,2) * t715;
t791 = t581 * t639;
t825 = (-m(3) * t791 + t575) * t617 / 0.2e1;
t666 = qJ(3,2) + t677;
t667 = -qJ(3,2) + t677;
t773 = sin(t666) + sin(t667);
t601 = t659 * t813 + sin(t689) * t818 + t656 * t817 - t773 * pkin(3);
t770 = cos(t666) + cos(t667);
t604 = t656 * t812 + cos(t689) * t818 + t659 * t817 - t770 * pkin(3);
t730 = 0.2e1 * qJ(3,2);
t699 = sin(t730);
t692 = pkin(3) * t699;
t714 = sin(qJ(3,2));
t616 = 0.1e1 / (t714 * t827 + t692 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1));
t589 = (t601 * t724 + t604 * t725) * t735 * t616;
t586 = t589 ^ 2;
t607 = t610 ^ 2;
t779 = t714 * t685;
t620 = 0.1e1 / (t692 + 0.2e1 * t779);
t702 = cos(t730);
t705 = t717 ^ 2;
t778 = t714 * t717;
t577 = (t607 * t717 * t815 + (-0.2e1 * t586 * t631 + ((-t631 * t702 + t705 * t814 - t685) * t610 + (t699 + 0.2e1 * t778) * t819 * t589) * t610) * pkin(3)) * t620;
t574 = t728 * t577;
t785 = t589 * t714;
t760 = pkin(3) * t785;
t784 = t589 * t717;
t580 = (-t819 * t760 + (t705 * t734 + t717 * t775 + t804 * t763 + t747) * t610) / (t692 / 0.2e1 + t779) * t610 * t735 + 0.2e1 * (t631 * t784 - t714 * t830) * t620 * t589;
t637 = rSges(3,1) * t717 - rSges(3,2) * t714;
t792 = t580 * t637;
t824 = (-m(3) * t792 + t574) * t616 / 0.2e1;
t662 = qJ(3,3) + t676;
t663 = -qJ(3,3) + t676;
t774 = sin(t662) + sin(t663);
t600 = t658 * t813 + sin(t688) * t818 + t655 * t817 - t774 * pkin(3);
t771 = cos(t662) + cos(t663);
t603 = t655 * t812 + cos(t688) * t818 + t658 * t817 - t771 * pkin(3);
t729 = 0.2e1 * qJ(3,3);
t698 = sin(t729);
t691 = pkin(3) * t698;
t713 = sin(qJ(3,3));
t615 = 0.1e1 / (t713 * t827 + t691 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1));
t588 = (t600 * t724 + t603 * t725) * t735 * t615;
t585 = t588 ^ 2;
t606 = t609 ^ 2;
t781 = t713 * t685;
t619 = 0.1e1 / (t691 + 0.2e1 * t781);
t701 = cos(t729);
t704 = t716 ^ 2;
t780 = t713 * t716;
t576 = (t606 * t716 * t815 + (-0.2e1 * t585 * t630 + ((-t630 * t701 + t704 * t814 - t685) * t609 + (t698 + 0.2e1 * t780) * t819 * t588) * t609) * pkin(3)) * t619;
t573 = t728 * t576;
t787 = t588 * t713;
t761 = pkin(3) * t787;
t786 = t588 * t716;
t579 = (-t819 * t761 + (t704 * t734 + t716 * t775 + t805 * t763 + t747) * t609) / (t691 / 0.2e1 + t781) * t609 * t735 + 0.2e1 * (t630 * t786 - t713 * t831) * t619 * t588;
t635 = rSges(3,1) * t716 - rSges(3,2) * t713;
t793 = t579 * t635;
t823 = (-m(3) * t793 + t573) * t615 / 0.2e1;
t582 = (-0.2e1 * t761 + t831) * t627 * t609;
t719 = rSges(3,3) + pkin(5);
t754 = t719 + t806;
t811 = m(3) * rSges(3,2);
t622 = -t754 * t811 + Icges(3,6);
t623 = t754 * rSges(3,1) * m(3) - Icges(3,5);
t612 = -t622 * t716 + t623 * t713;
t762 = rSges(3,1) * t811;
t686 = -Icges(3,4) + t762;
t722 = -Icges(3,1) / 0.2e1;
t732 = rSges(3,2) ^ 2;
t733 = rSges(3,1) ^ 2;
t767 = t732 + t733;
t738 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(3,2) / 0.2e1 + t722 - Icges(1,3) - Icges(2,3) + (-t737 + (0.2e1 * t806 - rSges(2,2)) * rSges(2,2) + (-0.2e1 * t694 - rSges(2,1)) * rSges(2,1)) * m(2) + (-t767 / 0.2e1 + (-0.2e1 * t806 - t719) * t719 + t828) * m(3);
t746 = t694 / 0.2e1 + pkin(2) / 0.2e1;
t741 = t609 * t746;
t745 = t806 / 0.4e1 + t719 / 0.4e1;
t766 = 0.2e1 * m(3);
t654 = (-t732 + t733) * m(3) + Icges(3,2) - Icges(3,1);
t807 = t654 / 0.2e1;
t808 = -t654 / 0.2e1;
t809 = Icges(3,6) / 0.4e1;
t810 = -Icges(3,5) / 0.4e1;
t822 = t627 * (0.4e1 * ((t713 * t809 + t716 * t810) * t588 + (t780 * t807 + (t704 - 0.1e1 / 0.2e1) * t686) * t609 + ((t716 * t741 - t745 * t787) * rSges(3,2) + (t713 * t741 + t745 * t786) * rSges(3,1)) * m(3)) * t588 - (t701 * t808 + t686 * t698 + (-t635 * t694 + (-t694 - t635) * pkin(2)) * t766 + t738) * t582 - t612 * t579);
t583 = (-0.2e1 * t760 + t830) * t628 * t610;
t613 = -t622 * t717 + t623 * t714;
t740 = t610 * t746;
t821 = t628 * (0.4e1 * ((t714 * t809 + t717 * t810) * t589 + (t778 * t807 + (t705 - 0.1e1 / 0.2e1) * t686) * t610 + ((t717 * t740 - t745 * t785) * rSges(3,2) + (t714 * t740 + t745 * t784) * rSges(3,1)) * m(3)) * t589 - (t702 * t808 + t686 * t699 + (-t637 * t694 + (-t694 - t637) * pkin(2)) * t766 + t738) * t583 - t613 * t580);
t584 = (-0.2e1 * t759 + t829) * t629 * t611;
t614 = -t622 * t718 + t623 * t715;
t739 = t746 * t611;
t820 = t629 * (0.4e1 * ((t715 * t809 + t718 * t810) * t590 + (t776 * t807 + (t706 - 0.1e1 / 0.2e1) * t686) * t611 + ((t718 * t739 - t745 * t783) * rSges(3,2) + (t715 * t739 + t745 * t782) * rSges(3,1)) * m(3)) * t590 - (t703 * t808 + t686 * t700 + (-t639 * t694 + (-t694 - t639) * pkin(2)) * t766 + t738) * t584 - t614 * t581);
t790 = t585 * (rSges(3,1) * t713 + rSges(3,2) * t716);
t789 = t586 * (rSges(3,1) * t714 + rSges(3,2) * t717);
t788 = t587 * (rSges(3,1) * t715 + rSges(3,2) * t718);
t768 = -t762 / 0.2e1 + Icges(3,4) / 0.2e1;
t765 = 0.2e1 * pkin(1);
t758 = t685 * t826;
t757 = t615 * t790;
t756 = t616 * t789;
t755 = t617 * t788;
t748 = rSges(3,1) * t758;
t626 = (t733 / 0.2e1 - t732 / 0.2e1) * m(3) + t722 + Icges(3,2) / 0.2e1;
t633 = rSges(3,2) * t758;
t673 = -t767 * m(3) - Icges(3,3);
t744 = (-m(3) * t576 * t635 + t579 * t673 + t582 * t612 + 0.2e1 * (t686 * t704 + (t626 * t713 + t633) * t716 + t713 * t748 + t768) * t606) * t615;
t743 = (-m(3) * t577 * t637 + t580 * t673 + t583 * t613 + 0.2e1 * (t686 * t705 + (t626 * t714 + t633) * t717 + t714 * t748 + t768) * t607) * t616;
t742 = (-m(3) * t578 * t639 + t581 * t673 + t584 * t614 + 0.2e1 * (t686 * t706 + (t626 * t715 + t633) * t718 + t715 * t748 + t768) * t608) * t617;
t684 = -qJ(3,1) + t690;
t683 = qJ(3,1) + t690;
t682 = -qJ(3,2) + t689;
t681 = qJ(3,2) + t689;
t680 = -qJ(3,3) + t688;
t679 = qJ(3,3) + t688;
t672 = -0.2e1 * qJ(3,1) + t678;
t669 = t731 + t678;
t668 = -0.2e1 * qJ(3,2) + t677;
t665 = t730 + t677;
t664 = -0.2e1 * qJ(3,3) + t676;
t661 = t729 + t676;
t596 = t769 * t827 + (cos(t684) + cos(t683)) * t765 + t772 * t813 + (cos(t672) + cos(t669) + 0.2e1 * t660) * pkin(3);
t595 = t770 * t827 + (cos(t682) + cos(t681)) * t765 + t773 * t813 + (cos(t668) + cos(t665) + 0.2e1 * t659) * pkin(3);
t594 = t771 * t827 + (cos(t680) + cos(t679)) * t765 + t774 * t813 + (cos(t664) + cos(t661) + 0.2e1 * t658) * pkin(3);
t593 = t772 * t827 + (sin(t684) + sin(t683)) * t765 + t769 * t812 + (sin(t672) + sin(t669) + 0.2e1 * t657) * pkin(3);
t592 = t773 * t827 + (sin(t682) + sin(t681)) * t765 + t770 * t812 + (sin(t668) + sin(t665) + 0.2e1 * t656) * pkin(3);
t591 = t774 * t827 + (sin(t680) + sin(t679)) * t765 + t771 * t812 + (sin(t664) + sin(t661) + 0.2e1 * t655) * pkin(3);
t1 = [t657 * t820 + t656 * t821 + t655 * t822 + (t603 * t744 + t604 * t743 + t605 * t742) * t735 + t594 * t823 + t595 * t824 + t596 * t825 + (-t594 * t757 - t595 * t756 - t596 * t755) * t826; -t660 * t820 - t659 * t821 - t658 * t822 + (t600 * t744 + t601 * t743 + t602 * t742) * t735 + t591 * t823 + t592 * t824 + t593 * t825 + (-t591 * t757 - t592 * t756 - t593 * t755) * t826; t573 + t574 + t575 + (-t788 - t789 - t790 - t791 - t792 - t793) * m(3);];
taucX  = t1;
