% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:10
% EndTime: 2020-08-06 18:28:10
% DurationCPUTime: 0.42s
% Computational Cost: add. (423->104), mult. (570->185), div. (36->7), fcn. (378->18), ass. (0->71)
t647 = m(2) + m(3);
t646 = pkin(1) + pkin(5);
t674 = -mrSges(1,2) + mrSges(2,3);
t631 = legFrame(3,2);
t621 = sin(t631);
t624 = cos(t631);
t607 = t624 * g(1) - t621 * g(2);
t635 = sin(qJ(1,3));
t641 = cos(qJ(1,3));
t601 = t641 * g(3) + t607 * t635;
t604 = t621 * g(1) + t624 * g(2);
t634 = sin(qJ(3,3));
t628 = 0.1e1 / t634;
t640 = cos(qJ(3,3));
t673 = ((-t601 * mrSges(3,1) + mrSges(3,2) * t604) * t640 + t634 * (mrSges(3,1) * t604 + t601 * mrSges(3,2))) * t628;
t632 = legFrame(2,2);
t622 = sin(t632);
t625 = cos(t632);
t608 = t625 * g(1) - t622 * g(2);
t637 = sin(qJ(1,2));
t643 = cos(qJ(1,2));
t602 = t643 * g(3) + t608 * t637;
t605 = t622 * g(1) + t625 * g(2);
t636 = sin(qJ(3,2));
t629 = 0.1e1 / t636;
t642 = cos(qJ(3,2));
t672 = ((-t602 * mrSges(3,1) + mrSges(3,2) * t605) * t642 + t636 * (mrSges(3,1) * t605 + t602 * mrSges(3,2))) * t629;
t633 = legFrame(1,2);
t623 = sin(t633);
t626 = cos(t633);
t609 = t626 * g(1) - t623 * g(2);
t639 = sin(qJ(1,1));
t645 = cos(qJ(1,1));
t603 = t645 * g(3) + t609 * t639;
t606 = t623 * g(1) + t626 * g(2);
t638 = sin(qJ(3,1));
t630 = 0.1e1 / t638;
t644 = cos(qJ(3,1));
t671 = ((-t603 * mrSges(3,1) + mrSges(3,2) * t606) * t644 + t638 * (mrSges(3,1) * t606 + t603 * mrSges(3,2))) * t630;
t618 = t634 * pkin(3) + qJ(2,3);
t615 = 0.1e1 / t618;
t670 = t601 * t615;
t619 = t636 * pkin(3) + qJ(2,2);
t616 = 0.1e1 / t619;
t669 = t602 * t616;
t620 = t638 * pkin(3) + qJ(2,1);
t617 = 0.1e1 / t620;
t668 = t603 * t617;
t649 = m(2) * pkin(1) + t646 * m(3) + mrSges(1,1) - mrSges(2,2) + mrSges(3,3);
t610 = t649 * g(3);
t652 = -t634 * mrSges(3,1) - mrSges(3,2) * t640 - t647 * qJ(2,3) - t674;
t667 = t615 * (t610 * t641 + (t635 * t649 + t652 * t641) * t607 - t635 * t652 * g(3));
t651 = -t636 * mrSges(3,1) - mrSges(3,2) * t642 - t647 * qJ(2,2) - t674;
t666 = t616 * (t610 * t643 + (t637 * t649 + t651 * t643) * t608 - t637 * t651 * g(3));
t650 = -t638 * mrSges(3,1) - mrSges(3,2) * t644 - t647 * qJ(2,1) - t674;
t665 = t617 * (t610 * t645 + (t639 * t649 + t650 * t645) * t609 - t639 * t650 * g(3));
t664 = t640 * qJ(2,3);
t663 = t642 * qJ(2,2);
t662 = t644 * qJ(2,1);
t661 = t628 * t670;
t660 = t629 * t669;
t659 = t630 * t668;
t658 = t641 * t667;
t657 = t643 * t666;
t656 = t645 * t665;
t627 = pkin(6) + t646;
t655 = qJ(2,1) * t639 + t627 * t645;
t654 = qJ(2,2) * t637 + t627 * t643;
t653 = qJ(2,3) * t635 + t627 * t641;
t648 = 0.1e1 / pkin(3);
t1 = [-g(1) * m(4) + t624 * t658 + t625 * t657 + t626 * t656 + (-t621 * t673 - t622 * t672 - t623 * t671) * t648 + (-(t655 * t626 * t638 + t623 * t662 + (t623 * t644 * t638 + (-t644 ^ 2 + 0.1e1) * t626 * t639) * pkin(3)) * t659 - (t654 * t625 * t636 + t622 * t663 + (t622 * t642 * t636 + (-t642 ^ 2 + 0.1e1) * t625 * t637) * pkin(3)) * t660 - (t653 * t624 * t634 + t621 * t664 + (t621 * t640 * t634 + (-t640 ^ 2 + 0.1e1) * t624 * t635) * pkin(3)) * t661) * t647; -t621 * t658 - t622 * t657 - t623 * t656 - g(2) * m(4) + (-t624 * t673 - t625 * t672 - t626 * t671) * t648 + (-((t626 * pkin(3) * t644 - t655 * t623) * t638 + t639 * pkin(3) * (t644 - 0.1e1) * (t644 + 0.1e1) * t623 + t626 * t662) * t659 - ((t625 * pkin(3) * t642 - t654 * t622) * t636 + t637 * pkin(3) * (t642 - 0.1e1) * (t642 + 0.1e1) * t622 + t625 * t663) * t660 - ((t624 * pkin(3) * t640 - t653 * t621) * t634 + t635 * pkin(3) * (t640 - 0.1e1) * (t640 + 0.1e1) * t621 + t624 * t664) * t661) * t647; -t635 * t667 - t637 * t666 - t639 * t665 - g(3) * m(4) + (-(t620 * t645 - t627 * t639) * t668 - (t619 * t643 - t627 * t637) * t669 - (t618 * t641 - t627 * t635) * t670) * t647;];
taugX  = t1;
