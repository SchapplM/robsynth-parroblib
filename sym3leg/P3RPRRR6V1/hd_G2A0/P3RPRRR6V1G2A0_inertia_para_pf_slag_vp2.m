% Calculate inertia matrix for parallel robot
% P3RPRRR6V1G2A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:17
% EndTime: 2020-08-06 18:36:19
% DurationCPUTime: 1.76s
% Computational Cost: add. (3774->255), mult. (4833->431), div. (531->10), fcn. (3690->53), ass. (0->173)
t657 = sin(pkin(7));
t677 = -pkin(6) - pkin(5);
t619 = t677 * t657 - pkin(1);
t664 = sin(qJ(1,3));
t604 = t619 * t664;
t658 = cos(pkin(7));
t670 = cos(qJ(1,3));
t717 = t670 * t677;
t718 = t670 * t657;
t753 = t604 - pkin(2) * t718 - (t664 * pkin(2) + t717) * t658;
t666 = sin(qJ(1,2));
t605 = t619 * t666;
t672 = cos(qJ(1,2));
t715 = t672 * t677;
t716 = t672 * t657;
t752 = t605 - pkin(2) * t716 - (t666 * pkin(2) + t715) * t658;
t668 = sin(qJ(1,1));
t606 = t619 * t668;
t674 = cos(qJ(1,1));
t713 = t674 * t677;
t714 = t674 * t657;
t751 = t606 - pkin(2) * t714 - (t668 * pkin(2) + t713) * t658;
t636 = t658 * pkin(1);
t627 = t636 + pkin(2);
t750 = -0.2e1 * pkin(1);
t749 = -0.2e1 * pkin(2);
t748 = 0.2e1 * pkin(2);
t747 = 0.2e1 * t677;
t746 = pkin(1) * t657;
t669 = cos(qJ(3,3));
t628 = t669 * pkin(3) + pkin(2);
t671 = cos(qJ(3,2));
t629 = t671 * pkin(3) + pkin(2);
t673 = cos(qJ(3,1));
t630 = t673 * pkin(3) + pkin(2);
t663 = sin(qJ(3,3));
t745 = t663 * mrSges(3,2);
t665 = sin(qJ(3,2));
t744 = t665 * mrSges(3,2);
t667 = sin(qJ(3,1));
t743 = t667 * mrSges(3,2);
t679 = 0.1e1 / pkin(3);
t727 = t628 * t664;
t742 = ((t717 + t727) * t658 - t604 + t628 * t718) * t679;
t726 = t629 * t666;
t741 = ((t715 + t726) * t658 - t605 + t629 * t716) * t679;
t725 = t630 * t668;
t740 = ((t713 + t725) * t658 - t606 + t630 * t714) * t679;
t647 = qJ(1,3) + pkin(7);
t633 = sin(t647);
t704 = pkin(7) + qJ(3,3);
t707 = -pkin(7) + qJ(3,3);
t739 = (t633 * t747 + cos(t647) * t749 + t670 * t750 + (-cos(qJ(1,3) - t707) - cos(qJ(1,3) + t704)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t663 * t748 + (sin(t704) + sin(t707)) * pkin(1));
t648 = qJ(1,2) + pkin(7);
t634 = sin(t648);
t705 = pkin(7) + qJ(3,2);
t708 = -pkin(7) + qJ(3,2);
t738 = (t634 * t747 + cos(t648) * t749 + t672 * t750 + (-cos(qJ(1,2) - t708) - cos(qJ(1,2) + t705)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t665 * t748 + (sin(t705) + sin(t708)) * pkin(1));
t649 = qJ(1,1) + pkin(7);
t635 = sin(t649);
t706 = pkin(7) + qJ(3,1);
t709 = -pkin(7) + qJ(3,1);
t737 = (t635 * t747 + cos(t649) * t749 + t674 * t750 + (-cos(qJ(1,1) - t709) - cos(qJ(1,1) + t706)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t667 * t748 + (sin(t706) + sin(t709)) * pkin(1));
t660 = legFrame(3,2);
t621 = t660 + t647;
t622 = -t660 + t647;
t736 = -sin(t621) / 0.2e1 + sin(t622) / 0.2e1;
t661 = legFrame(2,2);
t623 = t661 + t648;
t624 = -t661 + t648;
t735 = -sin(t623) / 0.2e1 + sin(t624) / 0.2e1;
t662 = legFrame(1,2);
t625 = t662 + t649;
t626 = -t662 + t649;
t734 = -sin(t625) / 0.2e1 + sin(t626) / 0.2e1;
t733 = cos(t622) / 0.2e1 + cos(t621) / 0.2e1;
t732 = cos(t624) / 0.2e1 + cos(t623) / 0.2e1;
t731 = cos(t626) / 0.2e1 + cos(t625) / 0.2e1;
t607 = 0.1e1 / (t636 + t628);
t651 = 0.1e1 / t663;
t730 = t607 * t651;
t608 = 0.1e1 / (t636 + t629);
t652 = 0.1e1 / t665;
t729 = t608 * t652;
t609 = 0.1e1 / (t636 + t630);
t653 = 0.1e1 / t667;
t728 = t609 * t653;
t637 = sin(t660);
t724 = t663 * t637;
t640 = cos(t660);
t723 = t663 * t640;
t638 = sin(t661);
t722 = t665 * t638;
t641 = cos(t661);
t721 = t665 * t641;
t639 = sin(t662);
t720 = t667 * t639;
t642 = cos(t662);
t719 = t667 * t642;
t712 = t627 * mrSges(3,1);
t711 = mrSges(3,2) * t749;
t710 = -0.2e1 * t636;
t654 = t669 ^ 2;
t703 = (t664 * t658 + t718) * t654 * pkin(3);
t655 = t671 ^ 2;
t702 = (t666 * t658 + t716) * t655 * pkin(3);
t656 = t673 ^ 2;
t701 = (t668 * t658 + t714) * t656 * pkin(3);
t700 = -m(3) * pkin(2) - mrSges(2,1);
t699 = t637 * t742;
t698 = t640 * t742;
t697 = ((t628 * t670 - t664 * t677) * t658 - t619 * t670 - t657 * t727) * t651 * t669;
t696 = t638 * t741;
t695 = t641 * t741;
t694 = ((t629 * t672 - t666 * t677) * t658 - t619 * t672 - t657 * t726) * t652 * t671;
t693 = t639 * t740;
t692 = t642 * t740;
t691 = ((t630 * t674 - t668 * t677) * t658 - t619 * t674 - t657 * t725) * t653 * t673;
t690 = t679 * t739;
t689 = t679 * t738;
t688 = t679 * t737;
t687 = pkin(5) + t746;
t614 = -t687 * mrSges(3,2) + Ifges(3,6);
t615 = t687 * mrSges(3,1) - Ifges(3,5);
t589 = t614 * t669 - t663 * t615;
t686 = t589 * t651 * t742;
t590 = t614 * t671 - t665 * t615;
t685 = t590 * t652 * t741;
t591 = t614 * t673 - t667 * t615;
t684 = t591 * t653 * t740;
t678 = m(2) + m(3);
t680 = Ifges(3,1) + Ifges(1,3) + Ifges(2,3) + 0.2e1 * (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t746 + (pkin(2) ^ 2 + pkin(5) ^ 2) * m(3) + t678 * pkin(1) ^ 2 + 0.2e1 * mrSges(3,3) * pkin(5);
t659 = -Ifges(3,1) + Ifges(3,2);
t618 = t673 * mrSges(3,1) - t743;
t617 = t671 * mrSges(3,1) - t744;
t616 = t669 * mrSges(3,1) - t745;
t579 = t659 * t656 + 0.2e1 * (Ifges(3,4) * t667 + t712) * t673 + (t700 + t743) * t710 + t667 * t711 + t680;
t578 = t659 * t655 + 0.2e1 * (Ifges(3,4) * t665 + t712) * t671 + (t700 + t744) * t710 + t665 * t711 + t680;
t577 = t659 * t654 + 0.2e1 * (Ifges(3,4) * t663 + t712) * t669 + (t700 + t745) * t710 + t663 * t711 + t680;
t576 = t642 * t701 + (pkin(3) * t720 - t751 * t642) * t673 + t627 * t720;
t575 = -t639 * t701 + (pkin(3) * t719 + t751 * t639) * t673 + t627 * t719;
t574 = t641 * t702 + (pkin(3) * t722 - t752 * t641) * t671 + t627 * t722;
t573 = -t638 * t702 + (pkin(3) * t721 + t752 * t638) * t671 + t627 * t721;
t572 = t640 * t703 + (pkin(3) * t724 - t753 * t640) * t669 + t627 * t724;
t571 = -t637 * t703 + (pkin(3) * t723 + t753 * t637) * t669 + t627 * t723;
t570 = t609 * t678 * t691 + t618 * t688;
t569 = t608 * t678 * t694 + t617 * t689;
t568 = t607 * t678 * t697 + t616 * t690;
t567 = -t635 * t609 * t579 + t591 * t688;
t566 = -t634 * t608 * t578 + t590 * t689;
t565 = -t633 * t607 * t577 + t589 * t690;
t564 = (t576 * t678 - t618 * t692) * t728;
t563 = (t575 * t678 + t618 * t693) * t728;
t562 = (t574 * t678 - t617 * t695) * t729;
t561 = (t573 * t678 + t617 * t696) * t729;
t560 = (t572 * t678 - t616 * t698) * t730;
t559 = (t571 * t678 + t616 * t699) * t730;
t558 = Ifges(3,3) * t688 + (-t591 * t635 + t618 * t691) * t609;
t557 = Ifges(3,3) * t689 + (-t590 * t634 + t617 * t694) * t608;
t556 = Ifges(3,3) * t690 + (-t589 * t633 + t616 * t697) * t607;
t555 = (t579 * t731 - t642 * t684) * t609;
t554 = (t578 * t732 - t641 * t685) * t608;
t553 = (t577 * t733 - t640 * t686) * t607;
t552 = (t579 * t734 + t639 * t684) * t609;
t551 = (t578 * t735 + t638 * t685) * t608;
t550 = (t577 * t736 + t637 * t686) * t607;
t549 = (t591 * t731 + (-Ifges(3,3) * t692 + t576 * t618) * t653) * t609;
t548 = (t590 * t732 + (-Ifges(3,3) * t695 + t574 * t617) * t652) * t608;
t547 = (t589 * t733 + (-Ifges(3,3) * t698 + t572 * t616) * t651) * t607;
t546 = (t591 * t734 + (Ifges(3,3) * t693 + t575 * t618) * t653) * t609;
t545 = (t590 * t735 + (Ifges(3,3) * t696 + t573 * t617) * t652) * t608;
t544 = (t589 * t736 + (Ifges(3,3) * t699 + t571 * t616) * t651) * t607;
t1 = [m(4) + (t555 * t731 + (-t549 * t692 + t564 * t576) * t653) * t609 + (t554 * t732 + (-t548 * t695 + t562 * t574) * t652) * t608 + (t553 * t733 + (-t547 * t698 + t560 * t572) * t651) * t607, (t555 * t734 + (t549 * t693 + t564 * t575) * t653) * t609 + (t554 * t735 + (t548 * t696 + t562 * t573) * t652) * t608 + (t553 * t736 + (t547 * t699 + t560 * t571) * t651) * t607, (-t555 * t635 + t564 * t691) * t609 + (-t554 * t634 + t562 * t694) * t608 + (-t553 * t633 + t560 * t697) * t607 + (t547 * t739 + t548 * t738 + t549 * t737) * t679; (t552 * t731 + (-t546 * t692 + t563 * t576) * t653) * t609 + (t551 * t732 + (-t545 * t695 + t561 * t574) * t652) * t608 + (t550 * t733 + (-t544 * t698 + t559 * t572) * t651) * t607, m(4) + (t552 * t734 + (t546 * t693 + t563 * t575) * t653) * t609 + (t551 * t735 + (t545 * t696 + t561 * t573) * t652) * t608 + (t550 * t736 + (t544 * t699 + t559 * t571) * t651) * t607, (-t552 * t635 + t563 * t691) * t609 + (-t551 * t634 + t561 * t694) * t608 + (-t550 * t633 + t559 * t697) * t607 + (t544 * t739 + t545 * t738 + t546 * t737) * t679; (t567 * t731 + (-t558 * t692 + t570 * t576) * t653) * t609 + (t566 * t732 + (-t557 * t695 + t569 * t574) * t652) * t608 + (t565 * t733 + (-t556 * t698 + t568 * t572) * t651) * t607, (t567 * t734 + (t558 * t693 + t570 * t575) * t653) * t609 + (t566 * t735 + (t557 * t696 + t569 * t573) * t652) * t608 + (t565 * t736 + (t556 * t699 + t568 * t571) * t651) * t607, m(4) + (-t567 * t635 + t570 * t691) * t609 + (-t566 * t634 + t569 * t694) * t608 + (-t565 * t633 + t568 * t697) * t607 + (t556 * t739 + t557 * t738 + t558 * t737) * t679;];
MX  = t1;
