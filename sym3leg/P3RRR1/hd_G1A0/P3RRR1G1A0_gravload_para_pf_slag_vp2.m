% Calculate Gravitation load for parallel robot
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRR1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:22
% EndTime: 2019-05-03 15:38:23
% DurationCPUTime: 0.31s
% Computational Cost: add. (459->85), mult. (614->153), div. (60->5), fcn. (536->20), ass. (0->75)
t669 = qJ(1,3) + qJ(2,3);
t654 = sin(t669);
t657 = cos(t669);
t675 = sin(qJ(1,3));
t678 = cos(qJ(1,3));
t700 = 0.1e1 / (t654 * t678 - t657 * t675);
t670 = qJ(1,2) + qJ(2,2);
t655 = sin(t670);
t658 = cos(t670);
t676 = sin(qJ(1,2));
t679 = cos(qJ(1,2));
t699 = 0.1e1 / (t655 * t679 - t658 * t676);
t671 = qJ(1,1) + qJ(2,1);
t656 = sin(t671);
t659 = cos(t671);
t677 = sin(qJ(1,1));
t680 = cos(qJ(1,1));
t698 = 0.1e1 / (t656 * t680 - t659 * t677);
t672 = legFrame(3,3);
t660 = sin(t672);
t663 = cos(t672);
t648 = -g(1) * t660 + g(2) * t663;
t651 = g(1) * t663 + g(2) * t660;
t612 = -t657 * (mrSges(2,1) * t648 - mrSges(2,2) * t651) + t654 * (mrSges(2,1) * t651 + mrSges(2,2) * t648);
t666 = m(2) * pkin(1) + mrSges(1,1);
t697 = ((mrSges(1,2) * t651 - t648 * t666) * t678 + t675 * (mrSges(1,2) * t648 + t651 * t666) + t612) * t700;
t673 = legFrame(2,3);
t661 = sin(t673);
t664 = cos(t673);
t649 = -g(1) * t661 + g(2) * t664;
t652 = g(1) * t664 + g(2) * t661;
t613 = -t658 * (mrSges(2,1) * t649 - mrSges(2,2) * t652) + t655 * (mrSges(2,1) * t652 + mrSges(2,2) * t649);
t696 = ((mrSges(1,2) * t652 - t649 * t666) * t679 + t676 * (mrSges(1,2) * t649 + t652 * t666) + t613) * t699;
t674 = legFrame(1,3);
t662 = sin(t674);
t665 = cos(t674);
t650 = -g(1) * t662 + g(2) * t665;
t653 = g(1) * t665 + g(2) * t662;
t614 = -t659 * (mrSges(2,1) * t650 - mrSges(2,2) * t653) + t656 * (mrSges(2,1) * t653 + mrSges(2,2) * t650);
t695 = ((mrSges(1,2) * t653 - t650 * t666) * t680 + t677 * (mrSges(1,2) * t650 + t653 * t666) + t614) * t698;
t694 = t612 * t700;
t693 = t613 * t699;
t692 = t614 * t698;
t627 = t654 * t663 + t657 * t660;
t628 = -t654 * t660 + t657 * t663;
t629 = t655 * t664 + t658 * t661;
t630 = -t655 * t661 + t658 * t664;
t631 = t656 * t665 + t659 * t662;
t632 = -t656 * t662 + t659 * t665;
t691 = 0.1e1 / pkin(1);
t690 = 0.1e1 / pkin(2);
t689 = koppelP(1,1);
t688 = koppelP(2,1);
t687 = koppelP(3,1);
t686 = koppelP(1,2);
t685 = koppelP(2,2);
t684 = koppelP(3,2);
t683 = mrSges(3,1);
t682 = mrSges(3,2);
t681 = xP(3);
t668 = cos(t681);
t667 = sin(t681);
t647 = -t667 * t686 + t668 * t689;
t646 = -t667 * t685 + t668 * t688;
t645 = -t667 * t684 + t668 * t687;
t644 = -t667 * t689 - t668 * t686;
t643 = -t667 * t688 - t668 * t685;
t642 = -t667 * t687 - t668 * t684;
t620 = pkin(1) * (-t662 * t677 + t665 * t680) + t632 * pkin(2);
t619 = pkin(1) * (-t661 * t676 + t664 * t679) + t630 * pkin(2);
t618 = pkin(1) * (-t660 * t675 + t663 * t678) + t628 * pkin(2);
t617 = pkin(1) * (t662 * t680 + t665 * t677) + t631 * pkin(2);
t616 = pkin(1) * (t661 * t679 + t664 * t676) + t629 * pkin(2);
t615 = pkin(1) * (t660 * t678 + t663 * t675) + t627 * pkin(2);
t1 = [-g(1) * m(3) + (t628 * t697 + t630 * t696 + t632 * t695 + (-t618 * t694 - t619 * t693 - t620 * t692) * t690) * t691; -g(2) * m(3) + (t627 * t697 + t629 * t696 + t631 * t695 + (-t615 * t694 - t616 * t693 - t617 * t692) * t690) * t691; -(-g(1) * t683 - g(2) * t682) * t667 + t668 * (g(1) * t682 - g(2) * t683) + ((t631 * t647 + t632 * t644) * t695 + (t629 * t646 + t630 * t643) * t696 + (t627 * t645 + t628 * t642) * t697 + (-(t617 * t647 + t620 * t644) * t692 - (t616 * t646 + t619 * t643) * t693 - (t615 * t645 + t618 * t642) * t694) * t690) * t691;];
taugX  = t1;
