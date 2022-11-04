% Return the minimum parameter vector for
% P6RRPRRR14V4G7A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,alpha3,alpha4,d1,d4,theta3]';
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
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 05:41
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P6RRPRRR14V4G7A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(6,3)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P6RRPRRR14V4G7A0_minimal_parameter_para: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V4G7A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V4G7A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6RRPRRR14V4G7A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V4G7A0_minimal_parameter_para: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
t12766 = sin(pkin(7));
t12768 = cos(pkin(7));
t12771 = Ifges(3,5) * t12766 + Ifges(3,6) * t12768;
t12767 = sin(pkin(3));
t12769 = cos(pkin(3));
t12777 = t12767 * t12769;
t12783 = t12771 * t12777;
t12763 = t12767 ^ 2;
t12765 = t12769 ^ 2;
t12762 = t12766 ^ 2;
t12764 = t12768 ^ 2;
t12774 = t12768 * t12766 * Ifges(3,4);
t12770 = Ifges(3,1) * t12762 + Ifges(3,2) * t12764 + 0.2e1 * t12774;
t12782 = t12763 * Ifges(3,3) + t12770 * t12765 + Ifges(2,2);
t12781 = Ifges(3,4) * (t12764 - t12762);
t12779 = Ifges(3,5) * t12768;
t12776 = 0.2e1 * t12783;
t12773 = (Ifges(3,1) - Ifges(3,2)) * t12768;
t1 = [Ifges(1,3) + t12782 - 0.2e1 * t12783; mrSges(1,1); mrSges(1,2) - mrSges(2,3); t12764 * Ifges(3,1) + t12762 * Ifges(3,2) + Ifges(2,1) - 0.2e1 * t12774 + t12776 - t12782; t12769 * t12781 - t12767 * t12779 + Ifges(2,4) + (Ifges(3,6) * t12767 + t12769 * t12773) * t12766; t12767 * t12781 + t12769 * t12779 + Ifges(2,5) + (-Ifges(3,6) * t12769 + t12767 * t12773) * t12766; Ifges(2,6) + t12771 * (t12765 - t12763) + (-Ifges(3,3) + t12770) * t12777; t12765 * Ifges(3,3) + t12770 * t12763 + Ifges(2,3) + t12776; mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,2); mrSges(3,3); m(3); Ifges(4,1); Ifges(4,4); Ifges(4,5); Ifges(4,2); Ifges(4,6); Ifges(4,3); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4);];
MPV  = t1;
