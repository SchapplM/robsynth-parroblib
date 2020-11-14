% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G3A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:04
% EndTime: 2020-08-06 19:10:08
% DurationCPUTime: 3.99s
% Computational Cost: add. (27090->275), mult. (41283->447), div. (6471->6), fcn. (34602->18), ass. (0->219)
t680 = cos(qJ(2,3));
t690 = pkin(1) + pkin(2);
t725 = t690 * t680;
t674 = sin(qJ(2,3));
t746 = t674 * qJ(3,3);
t640 = t725 + t746;
t675 = sin(qJ(1,3));
t681 = cos(qJ(1,3));
t616 = t681 * pkin(4) + t640 * t675;
t687 = xDP(3);
t804 = t616 * t687;
t682 = cos(qJ(2,2));
t724 = t690 * t682;
t676 = sin(qJ(2,2));
t743 = t676 * qJ(3,2);
t641 = t724 + t743;
t677 = sin(qJ(1,2));
t683 = cos(qJ(1,2));
t617 = t683 * pkin(4) + t641 * t677;
t803 = t617 * t687;
t684 = cos(qJ(2,1));
t723 = t690 * t684;
t678 = sin(qJ(2,1));
t740 = t678 * qJ(3,1);
t642 = t723 + t740;
t679 = sin(qJ(1,1));
t685 = cos(qJ(1,1));
t618 = t685 * pkin(4) + t642 * t679;
t802 = t618 * t687;
t778 = t675 * pkin(4);
t628 = t681 * t746 - t778;
t671 = legFrame(3,2);
t653 = sin(t671);
t656 = cos(t671);
t664 = t680 ^ 2;
t728 = t690 * t674;
t735 = t681 * t690;
t752 = t653 * qJ(3,3);
t595 = (t656 * t735 - t752) * t664 + (t628 * t656 + t653 * t728) * t680 + t752;
t749 = t656 * qJ(3,3);
t598 = (-t653 * t735 - t749) * t664 + (-t653 * t628 + t656 * t728) * t680 + t749;
t688 = xDP(2);
t689 = xDP(1);
t631 = 0.1e1 / t640;
t692 = 0.1e1 / qJ(3,3);
t755 = t631 * t692;
t580 = (t595 * t689 + t598 * t688 - t680 * t804) * t755;
t801 = t690 * t580;
t777 = t677 * pkin(4);
t629 = t683 * t743 - t777;
t672 = legFrame(2,2);
t654 = sin(t672);
t657 = cos(t672);
t665 = t682 ^ 2;
t727 = t690 * t676;
t732 = t683 * t690;
t751 = t654 * qJ(3,2);
t596 = (t657 * t732 - t751) * t665 + (t629 * t657 + t654 * t727) * t682 + t751;
t748 = t657 * qJ(3,2);
t599 = (-t654 * t732 - t748) * t665 + (-t654 * t629 + t657 * t727) * t682 + t748;
t632 = 0.1e1 / t641;
t694 = 0.1e1 / qJ(3,2);
t754 = t632 * t694;
t581 = (t596 * t689 + t599 * t688 - t682 * t803) * t754;
t800 = t690 * t581;
t776 = t679 * pkin(4);
t630 = t685 * t740 - t776;
t673 = legFrame(1,2);
t655 = sin(t673);
t658 = cos(t673);
t666 = t684 ^ 2;
t726 = t690 * t678;
t729 = t685 * t690;
t750 = t655 * qJ(3,1);
t597 = (t658 * t729 - t750) * t666 + (t630 * t658 + t655 * t726) * t684 + t750;
t747 = t658 * qJ(3,1);
t600 = (-t655 * t729 - t747) * t666 + (-t655 * t630 + t658 * t726) * t684 + t747;
t633 = 0.1e1 / t642;
t696 = 0.1e1 / qJ(3,1);
t753 = t633 * t696;
t582 = (t597 * t689 + t600 * t688 - t684 * t802) * t753;
t799 = t690 * t582;
t737 = t680 * qJ(3,3);
t701 = -t728 + t737;
t798 = t701 * t580;
t734 = t682 * qJ(3,2);
t702 = -t727 + t734;
t797 = t702 * t581;
t731 = t684 * qJ(3,1);
t703 = -t726 + t731;
t796 = t703 * t582;
t795 = t640 * t681 - t778;
t794 = t641 * t683 - t777;
t793 = t642 * t685 - t776;
t592 = (-t681 * t687 + (t653 * t688 - t656 * t689) * t675) * t631;
t691 = qJ(3,3) ^ 2;
t697 = pkin(4) ^ 2;
t607 = -t701 * t653 + t795 * t656;
t610 = -t795 * t653 - t701 * t656;
t583 = (t607 * t689 + t610 * t688 - t804) * t692;
t559 = -t583 + t801;
t745 = t674 * t559;
t553 = (-t580 * t737 + t745) * pkin(4) + ((qJ(3,3) + t690) * (-qJ(3,3) + t690) * t664 + 0.2e1 * t725 * t746 + t691 + t697) * t592;
t762 = t592 * t674;
t586 = pkin(4) * t762;
t565 = t586 + t801;
t736 = t680 * t692;
t763 = t592 * t664;
t764 = t592 * t631;
t541 = -t553 * t736 * t764 + ((-(t586 + t559) * t725 + (pkin(4) * t763 - t745) * qJ(3,3)) * t580 + (t565 * t680 + t580 * t746) * t583) * t755;
t699 = (pkin(1) ^ 2);
t710 = -t699 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t544 = (((-t691 + t710) * t580 + t690 * t583) * t580 + t565 * t583 + (pkin(4) * t798 - t553) * t592) * t692;
t785 = 0.2e1 * t583;
t550 = (-pkin(4) * t592 + t674 * t785 + 0.2e1 * t798) * t764;
t659 = m(3) * pkin(1) + mrSges(3,1);
t712 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t634 = t659 * qJ(3,3) + t712;
t568 = t634 * t580;
t577 = t580 ^ 2;
t771 = Ifges(2,6) - Ifges(3,6);
t646 = mrSges(3,2) * qJ(3,3) + t771;
t649 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t622 = -t646 * t680 + t674 * t649;
t772 = (-Ifges(2,1) - Ifges(3,1));
t779 = mrSges(3,1) * pkin(1);
t700 = m(3) * t699 + Ifges(2,2) + Ifges(3,3) + t772 + 2 * t779;
t786 = -2 * mrSges(3,3);
t714 = qJ(3,3) * t786;
t789 = -m(3) * t691 + t714;
t625 = t700 + t789;
t713 = -Ifges(1,3) + t772;
t744 = t674 * t680;
t650 = m(3) * qJ(3,3) + mrSges(3,3);
t767 = t583 * t650;
t770 = mrSges(3,2) * t674;
t782 = -0.2e1 * t634;
t722 = (-t625 * t664 + t744 * t782 + t713 + t789) * t550 + t622 * t541 - t544 * t770 + 0.4e1 * (t568 - t767 / 0.2e1) * t763 + (0.2e1 * (-t625 * t580 + t583 * t659) * t762 + t580 * (mrSges(3,2) * t785 - t649 * t580)) * t680 - t577 * t646 * t674 - 0.2e1 * (t568 - t767) * t592;
t792 = t675 * t722;
t593 = (-t683 * t687 + (t654 * t688 - t657 * t689) * t677) * t632;
t693 = qJ(3,2) ^ 2;
t608 = -t702 * t654 + t794 * t657;
t611 = -t794 * t654 - t702 * t657;
t584 = (t608 * t689 + t611 * t688 - t803) * t694;
t560 = -t584 + t800;
t742 = t676 * t560;
t554 = (-t581 * t734 + t742) * pkin(4) + ((qJ(3,2) + t690) * (-qJ(3,2) + t690) * t665 + 0.2e1 * t724 * t743 + t693 + t697) * t593;
t759 = t593 * t676;
t587 = pkin(4) * t759;
t566 = t587 + t800;
t733 = t682 * t694;
t760 = t593 * t665;
t761 = t593 * t632;
t542 = -t554 * t733 * t761 + ((-(t587 + t560) * t724 + (pkin(4) * t760 - t742) * qJ(3,2)) * t581 + (t566 * t682 + t581 * t743) * t584) * t754;
t545 = (((-t693 + t710) * t581 + t690 * t584) * t581 + t566 * t584 + (pkin(4) * t797 - t554) * t593) * t694;
t784 = 0.2e1 * t584;
t551 = (-pkin(4) * t593 + t676 * t784 + 0.2e1 * t797) * t761;
t635 = t659 * qJ(3,2) + t712;
t569 = t635 * t581;
t578 = t581 ^ 2;
t647 = mrSges(3,2) * qJ(3,2) + t771;
t623 = -t647 * t682 + t676 * t649;
t715 = qJ(3,2) * t786;
t788 = -m(3) * t693 + t715;
t626 = t700 + t788;
t741 = t676 * t682;
t651 = m(3) * qJ(3,2) + mrSges(3,3);
t766 = t584 * t651;
t769 = mrSges(3,2) * t676;
t781 = -0.2e1 * t635;
t721 = (-t626 * t665 + t741 * t781 + t713 + t788) * t551 + t623 * t542 - t545 * t769 + 0.4e1 * (t569 - t766 / 0.2e1) * t760 + (0.2e1 * (-t626 * t581 + t584 * t659) * t759 + t581 * (mrSges(3,2) * t784 - t649 * t581)) * t682 - t578 * t647 * t676 - 0.2e1 * (t569 - t766) * t593;
t791 = t677 * t721;
t594 = (-t685 * t687 + (t655 * t688 - t658 * t689) * t679) * t633;
t695 = qJ(3,1) ^ 2;
t609 = -t703 * t655 + t793 * t658;
t612 = -t793 * t655 - t703 * t658;
t585 = (t609 * t689 + t612 * t688 - t802) * t696;
t561 = -t585 + t799;
t739 = t678 * t561;
t555 = (-t582 * t731 + t739) * pkin(4) + ((qJ(3,1) + t690) * (-qJ(3,1) + t690) * t666 + 0.2e1 * t723 * t740 + t695 + t697) * t594;
t756 = t594 * t678;
t588 = pkin(4) * t756;
t567 = t588 + t799;
t730 = t684 * t696;
t757 = t594 * t666;
t758 = t594 * t633;
t543 = -t555 * t730 * t758 + ((-(t588 + t561) * t723 + (pkin(4) * t757 - t739) * qJ(3,1)) * t582 + (t567 * t684 + t582 * t740) * t585) * t753;
t546 = (((-t695 + t710) * t582 + t690 * t585) * t582 + t567 * t585 + (pkin(4) * t796 - t555) * t594) * t696;
t783 = 0.2e1 * t585;
t552 = (-pkin(4) * t594 + t678 * t783 + 0.2e1 * t796) * t758;
t636 = t659 * qJ(3,1) + t712;
t570 = t636 * t582;
t579 = t582 ^ 2;
t648 = qJ(3,1) * mrSges(3,2) + t771;
t624 = -t648 * t684 + t678 * t649;
t716 = qJ(3,1) * t786;
t787 = -m(3) * t695 + t716;
t627 = t700 + t787;
t738 = t678 * t684;
t652 = m(3) * qJ(3,1) + mrSges(3,3);
t765 = t585 * t652;
t768 = mrSges(3,2) * t678;
t780 = -0.2e1 * t636;
t720 = (-t627 * t666 + t738 * t780 + t713 + t787) * t552 + t624 * t543 - t546 * t768 + 0.4e1 * (t570 - t765 / 0.2e1) * t757 + (0.2e1 * (-t627 * t582 + t585 * t659) * t756 + t582 * (mrSges(3,2) * t783 - t649 * t582)) * t684 - t579 * t648 * t678 - 0.2e1 * (t570 - t765) * t594;
t790 = t679 * t720;
t589 = t592 ^ 2;
t711 = -2 * t779 - Ifges(3,2) - Ifges(2,3);
t719 = t622 * t550 + (-(t691 + t699) * m(3) + t714 + t711) * t541 + t659 * t544 + 0.2e1 * t580 * t767 + (t625 * t744 + t664 * t782 + t634) * t589;
t590 = t593 ^ 2;
t718 = t623 * t551 + (-(t693 + t699) * m(3) + t715 + t711) * t542 + t659 * t545 + 0.2e1 * t581 * t766 + (t626 * t741 + t665 * t781 + t635) * t590;
t591 = t594 ^ 2;
t717 = t624 * t552 + (-(t695 + t699) * m(3) + t716 + t711) * t543 + t659 * t546 + 0.2e1 * t582 * t765 + (t627 * t738 + t666 * t780 + t636) * t591;
t709 = t719 * t692;
t708 = t718 * t694;
t707 = t717 * t696;
t706 = (-m(3) * t544 + t659 * t541 - t550 * t770 - t577 * t650 + (t650 * t664 - t659 * t744 - t650) * t589) * t692;
t705 = (-m(3) * t545 + t659 * t542 - t551 * t769 - t578 * t651 + (t651 * t665 - t659 * t741 - t651) * t590) * t694;
t704 = (-m(3) * t546 + t659 * t543 - t552 * t768 - t579 * t652 + (t652 * t666 - t659 * t738 - t652) * t591) * t696;
t1 = [t609 * t704 + t608 * t705 + t607 * t706 + (t597 * t707 - t658 * t790) * t633 + (t596 * t708 - t657 * t791) * t632 + (t595 * t709 - t656 * t792) * t631; t612 * t704 + t611 * t705 + t610 * t706 + (t600 * t707 + t655 * t790) * t633 + (t599 * t708 + t654 * t791) * t632 + (t598 * t709 + t653 * t792) * t631; -t618 * t704 - t617 * t705 - t616 * t706 + (-t717 * t618 * t730 - t720 * t685) * t633 + (-t718 * t617 * t733 - t721 * t683) * t632 + (-t719 * t616 * t736 - t722 * t681) * t631;];
taucX  = t1;
